function Results41=run_rev41_monte_carlo(R)
%RUN_REV41_MONTE_CARLO Reproducible/resumable Rev.4.1 robustness ensemble.
%
% Default: 500 serial 120-s complete-6DOF trials.
%
% Smoke test:
%   R=rev41_mc_params(); R.N=20; R.resume=false;
%   Results41=run_rev41_monte_carlo(R);
%
% Final:
%   Results41=run_rev41_monte_carlo();

if nargin<1 || isempty(R), R=rev41_mc_params(); end

Pre=rev41_preflight(R);
Samples=rev41_generate_samples(R);

% Verify the theoretical/implementation initialization premise before any
% expensive trajectory integration.  Every uncertainty sample is retained.
InitialAudit=rev411_initial_admissibility_screen(Samples,R,Pre.Cnom);

outdir=fullfile(Pre.root,R.outputDir);
if ~exist(outdir,'dir'), mkdir(outdir); end
ckFile=fullfile(outdir,'REV41_CHECKPOINT.mat');

T=local_empty_run_template();
Runs=repmat(T,R.N,1);
startIdx=1;

if R.resume && exist(ckFile,'file')
    L=load(ckFile,'Runs','SamplesSaved','Rsave');
    same=(numel(L.Runs)==R.N) && (L.Rsave.seed==R.seed) && ...
         (numel(L.SamplesSaved)==R.N);
    if same
        Runs=L.Runs;
        Samples=L.SamplesSaved;
        done=find([Runs.index].'>0,1,'last');
        if ~isempty(done), startIdx=done+1; end
        fprintf('MC-RESUME checkpoint loaded; next trial = %d\n',startIdx);
    end
end

fprintf('\nREV.4.1 RANDOMIZED COMPLETE-6DOF ROBUSTNESS STUDY\n');
fprintf('=================================================\n');
fprintf('trials = %d, seed = %d\n',R.N,R.seed);
fprintf('serial resumable execution; controller gains remain frozen\n\n');

for k=startIdx:R.N
    q=rev41_run_one(Samples(k),R,Pre.Cnom);
    Runs(k)=q;

    if q.solverOK
        fprintf('MC %4d/%4d: hard=%d certDomain=%d latePos=%.3f m tube=%.3f deg lambdaMin=%.3f\n', ...
            k,R.N,q.hardFeasible,q.certifiedDomain,q.latePosNormMax,q.maxTubeDeg,q.minLambda);
    elseif q.domainExit
        fprintf('MC %4d/%4d: DYNAMIC DOMAIN EXIT: %s | %s\n', ...
            k,R.N,q.errorId,q.errorMessage);
    else
        fprintf('MC %4d/%4d: NUMERICAL/OTHER FAILURE: %s | %s\n', ...
            k,R.N,q.errorId,q.errorMessage);
    end

    if mod(k,R.checkpointEvery)==0 || k==R.N
        SamplesSaved=Samples; %#ok<NASGU>
        Rsave=R; %#ok<NASGU>
        save(ckFile,'Runs','SamplesSaved','Rsave','-v7.3');
    end
end

Sum=rev41_summarize(Runs,Samples,R);

Results41=struct();
Results41.status=Sum.status;
Results41.R=R;
Results41.Preflight=Pre;
Results41.InitialAudit=InitialAudit;
Results41.Samples=Samples;
Results41.Runs=Runs;
Results41.Summary=Sum;
Results41.Files=rev41_export_results(Results41);

fprintf('\nREV.4.1 MONTE-CARLO COMPLETION\n');
fprintf('==============================\n');
fprintf('full-run completion = %d/%d\n',Sum.Nsolver,R.N);
fprintf('dynamic governor/chart domain exits = %d/%d\n',Sum.NdomainExit,R.N);
fprintf('other numerical failures = %d/%d\n',Sum.NnumericalFailure,R.N);
fprintf('hard-feasibility preservation = %d/%d\n',Sum.Nhard,R.N);
fprintf('certified-domain preservation = %d/%d\n',Sum.Ncert,R.N);
fprintf('late position max quantiles [50 95 99]%% = [%.3f %.3f %.3f] m\n',Sum.posNormQuantiles);
fprintf('max attitude-tube quantiles [50 95 99]%% = [%.3f %.3f %.3f] deg\n',Sum.tubeQuantiles);
fprintf('FINAL STATUS: %s\n',Results41.status);

assignin('base','Results41',Results41);
end

function Q=local_empty_run_template()
% Must exactly match rev41_run_one's output field set.
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
