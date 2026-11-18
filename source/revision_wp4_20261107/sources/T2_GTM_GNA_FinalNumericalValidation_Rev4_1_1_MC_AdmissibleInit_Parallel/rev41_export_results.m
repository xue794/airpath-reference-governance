function Files=rev41_export_results(Results)
%REV41_EXPORT_RESULTS Write Monte-Carlo ledgers, summary, MAT, and paper plots.

root=fileparts(mfilename('fullpath'));
R=Results.R;
outdir=fullfile(root,R.outputDir);
if ~exist(outdir,'dir'), mkdir(outdir); end

summaryFile=fullfile(outdir,'REV41_MC_SUMMARY.txt');
samplesFile=fullfile(outdir,'REV41_MC_SAMPLES.csv');
runsFile=fullfile(outdir,'REV41_MC_RUNS.csv');
matFile=fullfile(outdir,'Rev41Results.mat');

local_write_summary(summaryFile,Results);
local_write_samples(samplesFile,Results.Samples);
local_write_runs(runsFile,Results.Runs);
save(matFile,'Results','-v7.3');

plots={};
if R.makePlots
    plots=local_plots(outdir,Results);
end

Files=struct('outputDir',outdir,'summary',summaryFile, ...
    'samples',samplesFile,'runs',runsFile,'mat',matFile,'plots',{plots});

fprintf('\nRev.4.1 exports written to:\n   %s\n',outdir);
end

function local_write_summary(file,Results)
fid=fopen(file,'w');
if fid<0, error('Rev41:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

R=Results.R; S=Results.Summary;
fprintf(fid,'REV.4.1 MONTE-CARLO ROBUSTNESS SUMMARY\n');
fprintf(fid,'======================================\n\n');
fprintf(fid,'STATUS: %s\n',S.status);
fprintf(fid,'Trials / seed: %d / %d\n',R.N,R.seed);
fprintf(fid,'Full-run completion: %d/%d (%.4f)\n',S.Nsolver,R.N,S.solverRate);
fprintf(fid,'Dynamic governor/chart domain exits: %d/%d (%.4f)\n',S.NdomainExit,R.N,S.domainExitRate);
fprintf(fid,'Other numerical failures: %d/%d (%.4f)\n',S.NnumericalFailure,R.N,S.numericalFailureRate);
fprintf(fid,'Hard-feasibility preservation: %d/%d (%.4f)\n',S.Nhard,R.N,S.hardFeasibleRate);
fprintf(fid,'Certified-domain preservation: %d/%d (%.4f)\n\n',S.Ncert,R.N,S.certifiedDomainRate);

fprintf(fid,'Late position-norm max quantiles [50 95 99] %%: %.9g %.9g %.9g m\n',S.posNormQuantiles);
fprintf(fid,'Attitude-tube max quantiles [50 95 99] %%: %.9g %.9g %.9g deg\n',S.tubeQuantiles);
fprintf(fid,'Minimum lambda quantiles [1 5 50] %%: %.9g %.9g %.9g\n',S.minLambdaQuantiles);
fprintf(fid,'Minimum phi quantiles [1 5 50] %%: %.9g %.9g %.9g\n',S.minPhiQuantiles);
fprintf(fid,'Late wind-error max quantiles [50 95 99] %%: %.9g %.9g %.9g m/s\n',S.windErrorQuantiles);
fprintf(fid,'Max physical moment-residual quantiles [50 95 99] %%: %.9g %.9g %.9g N m\n\n',S.dTauQuantiles);

fprintf(fid,'Worst late-position trial: %d\n',S.worstPositionIndex);
fprintf(fid,'Worst attitude-tube trial: %d\n',S.worstTubeIndex);
fprintf(fid,'Minimum-lambda trial: %d\n',S.minLambdaIndex);

if isfinite(S.pUpper95IfZeroHardFailures)
    fprintf(fid,'Conditional 95%% zero-failure upper probability under THIS STUDY SAMPLING MODEL: %.9g\n', ...
        S.pUpper95IfZeroHardFailures);
end

fprintf(fid,'\nINTERPRETATION LOCK\n');
fprintf(fid,'-------------------\n');
fprintf(fid,'1. Rev.4.1 is an extended empirical robustness study, not a stability proof.\n');
fprintf(fid,'2. The uncertainty distributions are study-defined, not identified physical probability laws.\n');
fprintf(fid,'3. Controller/governor/allocator gains remain frozen; only plant/environment/actuator/initial conditions vary.\n');
fprintf(fid,'4. Failure rates, if reported, are conditional on this artificial sampling ensemble and are not real-world reliability estimates.\n');
end

function local_write_samples(file,S)
fid=fopen(file,'w');
if fid<0, error('Rev41:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

fprintf(fid,['index,massScale,IxxScale,IyyScale,IzzScale,IxzScale,' ...
    'servo1TauScale,servo2TauScale,servo3TauScale,thrustTauScale,' ...
    'windMeanN,windMeanE,windMeanD,windAmpN,windAmpE,windAmpD,' ...
    'windOmegaN,windOmegaE,windOmegaD,noiseScale,' ...
    'initPosN,initPosE,initPosD,initVgBu,initVgBv,initVgBw,' ...
    'initAttRotNormDeg,initOmegaPdegps,initOmegaQdegps,initOmegaRdegps,' ...
    'initWhatN,initWhatE,initWhatD,windPeakBound,aeroMaxDeviationPct\n']);

for k=1:numel(S)
    q=S(k);
    v=[q.index q.massScale q.inertiaScale.' q.servoTauScale.' q.thrustTauScale ...
       q.windMean.' q.windAmp.' q.windOmega.' q.noiseScale ...
       q.initPosition.' q.initVgB.' rad2deg(norm(q.initAttitudeRotvec)) ...
       rad2deg(q.initOmega).' q.initWhatError.' q.windPeakBound q.aeroMaxDeviationPct];
    fprintf(fid,'%.16e',v(1));
    fprintf(fid,',%.16e',v(2:end));
    fprintf(fid,'\n');
end
end

function local_write_runs(file,Runs)
fid=fopen(file,'w');
if fid<0, error('Rev41:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

fprintf(fid,['index,solverOK,hardFeasible,certifiedDomain,domainExit,numericalFailure,' ...
    'lateErN,lateErE,lateErD,latePosNormMax,latePosNormRMS,maxPosNorm,' ...
    'lateWindErrorMax,minLambda,minPhi,maxTubeDeg,Vmin,Vmax,gammaAbsMaxDeg,' ...
    'alphaMinDeg,alphaMaxDeg,betaAbsMaxDeg,kappaMin,Tmin,Tmax,maxThrustRate,' ...
    'delta1MaxDeg,delta2MaxDeg,delta3MaxDeg,' ...
    'deltaRate1MaxDegS,deltaRate2MaxDegS,deltaRate3MaxDegS,' ...
    'allocatorActiveSamples,maxDTau,lateDTau,maxSurfaceForce,' ...
    'maxSlideIdentity,maxWindKinematic,maxQNormError,maxQrNormError\n']);

for k=1:numel(Runs)
    q=Runs(k);
    v=[q.index q.solverOK q.hardFeasible q.certifiedDomain q.domainExit q.numericalFailure ...
       q.lateErLinf q.latePosNormMax q.latePosNormRMS q.maxPosNorm ...
       q.lateWindErrorMax q.minLambda q.minPhi q.maxTubeDeg q.Vmin q.Vmax ...
       q.gammaAbsMaxDeg q.alphaMinDeg q.alphaMaxDeg q.betaAbsMaxDeg q.kappaMin ...
       q.Tmin q.Tmax q.maxThrustRate q.maxDeltaActDeg q.maxDeltaRateDegS ...
       q.allocatorActiveSamples q.maxDTau q.lateDTau q.maxSurfaceForce ...
       q.maxSlideIdentity q.maxWindKinematic q.maxQNormError q.maxQrNormError];
    fprintf(fid,'%.16e',v(1));
    fprintf(fid,',%.16e',v(2:end));
    fprintf(fid,'\n');
end
end

function P=local_plots(outdir,Results)
Runs=Results.Runs;
Samples=Results.Samples;
R=Results.R;

ok=[Runs.solverOK].';
idx=find(ok);
if isempty(idx), P={}; return; end

pos=[Runs(idx).latePosNormMax].';
tube=[Runs(idx).maxTubeDeg].';
lam=[Runs(idx).minLambda].';
wind=[Samples(idx).windPeakBound].';
aero=[Samples(idx).aeroMaxDeviationPct].';

P={};

% Fig. 1: empirical CDFs.
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.4]);
subplot(2,1,1);
[x,y]=local_ecdf(pos); plot(x,y,'LineWidth',1.2); grid on; box on;
xlabel('late max ||e_r|| [m]','Interpreter','none'); ylabel('empirical CDF');
subplot(2,1,2);
[x,y]=local_ecdf(tube); plot(x,y,'LineWidth',1.2); grid on; box on;
xlabel('max attitude tube [deg]','Interpreter','none'); ylabel('empirical CDF');
P{end+1}=local_export(f,outdir,'Fig_Rev41_MC_CDF'); close(f); %#ok<AGROW>

% Fig. 2: governor distribution.
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 4.2]);
[x,y]=local_ecdf(lam); plot(x,y,'LineWidth',1.2); grid on; box on;
xlabel('minimum lambda*','Interpreter','none'); ylabel('empirical CDF');
P{end+1}=local_export(f,outdir,'Fig_Rev41_MC_GovernorCDF'); close(f); %#ok<AGROW>

% Fig. 3: severity scatter.
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.4]);
subplot(2,1,1);
plot(wind,pos,'o','MarkerSize',4); grid on; box on;
xlabel('wind componentwise peak-norm bound [m/s]','Interpreter','none');
ylabel('late max ||e_r|| [m]','Interpreter','none');
subplot(2,1,2);
plot(aero,pos,'o','MarkerSize',4); grid on; box on;
xlabel('max aero coefficient deviation [%]','Interpreter','none');
ylabel('late max ||e_r|| [m]','Interpreter','none');
P{end+1}=local_export(f,outdir,'Fig_Rev41_MC_SeverityScatter'); close(f); %#ok<AGROW>

% Fig. 4: worst-position trial trajectory, rerun once.
iw=Results.Summary.worstPositionIndex;
C=rev41_apply_sample(Results.Preflight.Cnom,Samples(iw),R);
S=simulate_complete_6dof_case(R.mode,C);

f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.6]);
labs={'e_N [m]','e_E [m]','e_D [m]'};
for j=1:3
    subplot(3,1,j);
    plot(S.t,S.Y.er(:,j),'LineWidth',1.15); grid on; box on;
    ylabel(labs{j},'Interpreter','none');
    if j==3, xlabel('t [s]','Interpreter','none'); end
end
P{end+1}=local_export(f,outdir,'Fig_Rev41_MC_WorstPositionTrial'); close(f); %#ok<AGROW>
end

function [x,y]=local_ecdf(v)
x=sort(v(isfinite(v)));
n=numel(x);
if n==0, y=[]; return; end
y=(1:n).'/n;
end

function F=local_export(fig,outdir,stem)
png=fullfile(outdir,[stem '.png']);
epsf=fullfile(outdir,[stem '.eps']);
print(fig,png,'-dpng','-r300');
print(fig,epsf,'-depsc2','-r600');
F=struct('png',png,'eps',epsf);
end
