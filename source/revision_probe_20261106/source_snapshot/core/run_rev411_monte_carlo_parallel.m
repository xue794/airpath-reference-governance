function Results411=run_rev411_monte_carlo_parallel(numWorkers)
%RUN_REV411_MONTE_CARLO_PARALLEL
% Parallel/resumable Rev.4.1.1 Monte-Carlo robustness study for MATLAB R2018a.
%
% Usage:
%   Results411 = run_rev411_monte_carlo_parallel();     % default 48 workers
%   Results411 = run_rev411_monte_carlo_parallel(56);   % use 56 workers
%
% The random ensemble is generated ONCE on the MATLAB client using the same
% Rev.4.1.1 seed and generator. No random numbers are drawn inside PARFOR.
% Therefore each trial keeps exactly the same frozen sample index and
% uncertainty realization as the serial runner.
%
% Checkpointing:
%   Results are committed after each parallel batch. If MATLAB is interrupted,
%   rerun the same command and completed batches are reused.
%
% No theory/controller/plant uncertainty range is changed by this file.

if nargin<1 || isempty(numWorkers)
    numWorkers=48;
end

rev411_package_guard();

R=rev41_mc_params();
Pre=rev41_preflight(R);

% Same frozen ensemble as the serial Rev.4.1.1 runner.
Samples=rev41_generate_samples(R);

% Mandatory assumption-consistency screen before expensive integration.
InitialAudit=rev411_initial_admissibility_screen(Samples,R,Pre.Cnom);

% Parallel output is kept separate from any serial Rev.4.1.1 results.
R.outputDir=sprintf('rev411_parallel_outputs_%dworkers',numWorkers);
outdir=fullfile(Pre.root,R.outputDir);
if ~exist(outdir,'dir'), mkdir(outdir); end
ckFile=fullfile(outdir,'REV411_PARALLEL_CHECKPOINT.mat');

% Open/reconfigure local process pool.
pool=gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers~=numWorkers
    delete(pool);
    pool=[];
end
if isempty(pool)
    fprintf('\nOpening local parallel pool with %d workers...\n',numWorkers);
    pool=parpool('local',numWorkers);
end

actualWorkers=pool.NumWorkers;
if actualWorkers~=numWorkers
    warning('Rev411Parallel:Workers', ...
        'Requested %d workers but pool has %d.',numWorkers,actualWorkers);
end

% Two worker-waves per checkpoint batch keeps workers fed despite unequal
% trial durations while limiting lost work after interruption.
batchSize=max(2*actualWorkers,actualWorkers+1);

T=local_empty_run_template();
Runs=repmat(T,R.N,1);

% Resume only if the frozen study identity matches.
if R.resume && exist(ckFile,'file')
    L=load(ckFile,'RunsSaved','SamplesSaved','Rsave','numWorkersSaved');
    same=numel(L.RunsSaved)==R.N && ...
         numel(L.SamplesSaved)==R.N && ...
         L.Rsave.seed==R.seed && ...
         L.Rsave.N==R.N && ...
         L.numWorkersSaved==actualWorkers;
    if same
        Runs=L.RunsSaved;
        Samples=L.SamplesSaved;
        fprintf('Parallel checkpoint loaded.\n');
    else
        warning('Rev411Parallel:Checkpoint', ...
            'Existing checkpoint identity differs; starting a fresh run.');
    end
end

done=[Runs.index].'>0;
pending=find(~done);

fprintf('\nREV.4.1.1 PARALLEL RANDOMIZED COMPLETE-6DOF STUDY\n');
fprintf('=================================================\n');
fprintf('trials / seed = %d / %d\n',R.N,R.seed);
fprintf('workers = %d\n',actualWorkers);
fprintf('parallel checkpoint batch size = %d trials\n',batchSize);
fprintf('already completed = %d, pending = %d\n',sum(done),numel(pending));
fprintf('controller gains and uncertainty draws remain frozen\n\n');

wallStart=tic;

for ib=1:batchSize:numel(pending)
    idx=pending(ib:min(ib+batchSize-1,numel(pending)));
    nBatch=numel(idx);

    fprintf('PARALLEL BATCH: trials %d...%d (%d trials)\n', ...
        idx(1),idx(end),nBatch);

    % Preallocate the exact output structure. Each PARFOR iteration writes
    % one independent sliced element.
    BatchRuns=repmat(T,nBatch,1);

    batchStart=tic;
    parfor jj=1:nBatch
        k=idx(jj);
        BatchRuns(jj)=rev41_run_one(Samples(k),R,Pre.Cnom);
    end
    batchWall=toc(batchStart);

    % Commit results serially, then checkpoint.
    for jj=1:nBatch
        k=idx(jj);
        Runs(k)=BatchRuns(jj);

        q=BatchRuns(jj);
        if q.solverOK
            fprintf('MC %4d/%4d: hard=%d certDomain=%d latePos=%.3f m tube=%.3f deg lambdaMin=%.3f\n', ...
                k,R.N,q.hardFeasible,q.certifiedDomain,q.latePosNormMax, ...
                q.maxTubeDeg,q.minLambda);
        elseif q.domainExit
            fprintf('MC %4d/%4d: DYNAMIC DOMAIN EXIT: %s | %s\n', ...
                k,R.N,q.errorId,q.errorMessage);
        else
            fprintf('MC %4d/%4d: NUMERICAL/OTHER FAILURE: %s | %s\n', ...
                k,R.N,q.errorId,q.errorMessage);
        end
    end

    RunsSaved=Runs; %#ok<NASGU>
    SamplesSaved=Samples; %#ok<NASGU>
    Rsave=R; %#ok<NASGU>
    numWorkersSaved=actualWorkers; %#ok<NASGU>
    save(ckFile,'RunsSaved','SamplesSaved','Rsave','numWorkersSaved','-v7.3');

    fprintf('Batch wall time = %.1f s; total completed = %d/%d\n\n', ...
        batchWall,sum([Runs.index].'>0),R.N);
end

totalWall=toc(wallStart);

Sum=rev41_summarize(Runs,Samples,R);

Results411=struct();
Results411.status=Sum.status;
Results411.R=R;
Results411.Preflight=Pre;
Results411.InitialAudit=InitialAudit;
Results411.Samples=Samples;
Results411.Runs=Runs;
Results411.Summary=Sum;
Results411.Parallel=struct( ...
    'numWorkers',actualWorkers, ...
    'batchSize',batchSize, ...
    'wallTimeSeconds',totalWall, ...
    'checkpointFile',ckFile);
Results411.Files=rev41_export_results(Results411);

fprintf('\nREV.4.1.1 PARALLEL MONTE-CARLO COMPLETION\n');
fprintf('==========================================\n');
fprintf('workers = %d\n',actualWorkers);
fprintf('parallel wall time = %.1f s (%.2f min)\n',totalWall,totalWall/60);
fprintf('full-run completion = %d/%d\n',Sum.Nsolver,R.N);
fprintf('dynamic governor/chart domain exits = %d/%d\n',Sum.NdomainExit,R.N);
fprintf('other numerical failures = %d/%d\n',Sum.NnumericalFailure,R.N);
fprintf('hard-feasibility preservation = %d/%d\n',Sum.Nhard,R.N);
fprintf('certified-domain preservation = %d/%d\n',Sum.Ncert,R.N);
fprintf('FINAL STATUS: %s\n',Results411.status);

assignin('base','Results411',Results411);
end

function Q=local_empty_run_template()
Q=struct( ...
    'index',0,'solverOK',false,'hardFeasible',false,'certifiedDomain',false, ...
    'domainExit',false,'numericalFailure',false, ...
    'lateErLinf',[nan nan nan],'latePosNormMax',nan,'latePosNormRMS',nan, ...
    'maxPosNorm',nan,'lateWindErrorMax',nan,'minLambda',nan,'minPhi',nan, ...
    'maxTubeDeg',nan,'Vmin',nan,'Vmax',nan,'gammaAbsMaxDeg',nan, ...
    'alphaMinDeg',nan,'alphaMaxDeg',nan,'betaAbsMaxDeg',nan,'kappaMin',nan, ...
    'Tmin',nan,'Tmax',nan,'maxThrustRate',nan, ...
    'maxDeltaActDeg',[nan nan nan],'maxDeltaRateDegS',[nan nan nan], ...
    'allocatorActiveSamples',nan,'maxDTau',nan,'lateDTau',nan, ...
    'maxSurfaceForce',nan,'maxSlideIdentity',nan,'maxWindKinematic',nan, ...
    'maxQNormError',nan,'maxQrNormError',nan, ...
    'flagGovernor',false,'flagAirspeed',false,'flagGamma',false, ...
    'flagAlpha',false,'flagSurfaceBox',false,'flagSurfaceRate',false, ...
    'flagThrustBox',false,'flagThrustRate',false,'flagQuaternion',false, ...
    'flagKappaPositive',false,'flagAttitudeTube',false,'flagBetaTube',false, ...
    'flagKappaCert',false,'errorId','','errorMessage','');
end
