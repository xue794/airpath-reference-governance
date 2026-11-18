function B=benchmark_rev411_parallel(workerList,nTrials)
%BENCHMARK_REV411_PARALLEL Quick scaling benchmark for Rev.4.1.1.
%
% Example:
%   B=benchmark_rev411_parallel([16 32 48 56],56);
%
% Uses the same frozen uncertainty samples and only the first nTrials.
% Results are timing diagnostics only; do not use this benchmark as the
% manuscript Monte-Carlo result.

if nargin<1 || isempty(workerList)
    workerList=[16 32 48 56];
end
if nargin<2 || isempty(nTrials)
    nTrials=56;
end

R=rev41_mc_params();
Pre=rev41_preflight(R);
Samples=rev41_generate_samples(R);
rev411_initial_admissibility_screen(Samples,R,Pre.Cnom);

nTrials=min(nTrials,R.N);
B=struct('workers',num2cell(workerList(:)), ...
         'wallSeconds',num2cell(nan(numel(workerList),1)), ...
         'trialsPerSecond',num2cell(nan(numel(workerList),1)));

for iw=1:numel(workerList)
    nw=workerList(iw);

    pool=gcp('nocreate');
    if ~isempty(pool), delete(pool); end
    pool=parpool('local',nw);

    T=local_empty_run_template();
    Q=repmat(T,nTrials,1);

    fprintf('\nBenchmark: %d workers, %d trials\n',nw,nTrials);
    tic;
    parfor k=1:nTrials
        Q(k)=rev41_run_one(Samples(k),R,Pre.Cnom);
    end
    sec=toc;

    B(iw).wallSeconds=sec;
    B(iw).trialsPerSecond=nTrials/sec;

    fprintf('wall = %.1f s, throughput = %.4f trials/s\n', ...
        sec,B(iw).trialsPerSecond);

    delete(pool);
end
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
