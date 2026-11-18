function S=rev41_summarize(Runs,Samples,R)
%REV41_SUMMARIZE Aggregate Rev.4.1 randomized robustness results.

ok=[Runs.solverOK].';
hard=[Runs.hardFeasible].';
cert=[Runs.certifiedDomain].';
domainExit=[Runs.domainExit].';
numericalFailure=[Runs.numericalFailure].';

S.N=R.N;
S.Nsolver=sum(ok);
S.Nhard=sum(hard);
S.Ncert=sum(cert);
S.NdomainExit=sum(domainExit);
S.NnumericalFailure=sum(numericalFailure);
S.solverRate=S.Nsolver/R.N;
S.hardFeasibleRate=S.Nhard/R.N;
S.certifiedDomainRate=S.Ncert/R.N;
S.domainExitRate=S.NdomainExit/R.N;
S.numericalFailureRate=S.NnumericalFailure/R.N;

if S.Nhard>0
    idx=find(hard);
else
    idx=find(ok);
end

if isempty(idx)
    S.status='REV41_NO_SUCCESSFUL_RUNS';
    S.pUpper95IfZeroHardFailures=nan;
    return;
end

pos=[Runs(idx).latePosNormMax].';
tube=[Runs(idx).maxTubeDeg].';
lam=[Runs(idx).minLambda].';
phi=[Runs(idx).minPhi].';
wd=[Runs(idx).lateWindErrorMax].';
dtau=[Runs(idx).maxDTau].';

S.posNormQuantiles=local_q(pos,[0.50 0.95 0.99]);
S.tubeQuantiles=local_q(tube,[0.50 0.95 0.99]);
S.minLambdaQuantiles=local_q(lam,[0.01 0.05 0.50]);
S.minPhiQuantiles=local_q(phi,[0.01 0.05 0.50]);
S.windErrorQuantiles=local_q(wd,[0.50 0.95 0.99]);
S.dTauQuantiles=local_q(dtau,[0.50 0.95 0.99]);

% Worst-case indices are selected only from successfully completed runs so
% solver-failure NaNs cannot contaminate the ranking.
okIdx=find(ok);
v=[Runs(okIdx).latePosNormMax].';
[~,jj]=max(v);
S.worstPositionIndex=Runs(okIdx(jj)).index;

v=[Runs(okIdx).maxTubeDeg].';
[~,jj]=max(v);
S.worstTubeIndex=Runs(okIdx(jj)).index;

v=[Runs(okIdx).minLambda].';
[~,jj]=min(v);
S.minLambdaIndex=Runs(okIdx(jj)).index;

nFail=R.N-S.Nhard;
if nFail==0
    S.pUpper95IfZeroHardFailures=1-0.05^(1/R.N);
else
    S.pUpper95IfZeroHardFailures=nan;
end

if S.NnumericalFailure>0
    S.status='REV411_MONTE_CARLO_COMPLETED_WITH_NUMERICAL_FAILURES';
elseif S.NdomainExit>0
    S.status='REV411_MONTE_CARLO_COMPLETED_WITH_DYNAMIC_DOMAIN_EXITS';
elseif S.Nsolver==R.N && S.Nhard==R.N
    S.status='REV411_MONTE_CARLO_ALL_TRIALS_COMPLETED_HARD_FEASIBLE';
else
    S.status='REV411_MONTE_CARLO_COMPLETED_WITH_FEASIBILITY_EXCEEDANCES';
end

% Study-ensemble descriptors.
S.windPeakBoundRange=[min([Samples.windPeakBound]) max([Samples.windPeakBound])];
S.aeroMaxDeviationPctRange=[min([Samples.aeroMaxDeviationPct]) max([Samples.aeroMaxDeviationPct])];
end

function y=local_q(x,p)
x=x(isfinite(x));
x=sort(x(:));
n=numel(x);
y=nan(size(p));
if n==0, return; end
for k=1:numel(p)
    h=1+(n-1)*p(k);
    lo=floor(h); hi=ceil(h);
    if lo==hi
        y(k)=x(lo);
    else
        y(k)=x(lo)+(h-lo)*(x(hi)-x(lo));
    end
end
end
