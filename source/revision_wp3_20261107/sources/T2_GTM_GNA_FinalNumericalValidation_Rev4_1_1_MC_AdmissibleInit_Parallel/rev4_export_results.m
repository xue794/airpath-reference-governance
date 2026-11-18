function Files=rev4_export_results(Results)
%REV4_EXPORT_RESULTS Export numerical summary, CSV ledgers, MAT file, figures.

root=fileparts(mfilename('fullpath'));
R=Results.Preflight.R;
outdir=fullfile(root,R.outputDir);

if ~exist(outdir,'dir'), mkdir(outdir); end

summaryFile=fullfile(outdir,'REV4_RESULTS_SUMMARY.txt');
caseFile=fullfile(outdir,'REV4_CASE_METRICS.csv');
macroFile=fullfile(outdir,'REV4_MACRO_CERTIFICATE.csv');
pFile=fullfile(outdir,'REV4_MACRO_STORAGE_P_NORMALIZED.csv');
matFile=fullfile(outdir,'Rev4Results.mat');

local_write_summary(summaryFile,Results);
local_write_cases(caseFile,Results.Extended.cases);
local_write_macro(macroFile,Results.Certificate.Family);
dlmwrite(pFile,Results.Certificate.P,'precision','%.16e');

if R.saveMat
    save(matFile,'Results','-v7.3');
end

plotFiles={};
if R.makePlots
    plotFiles=local_make_plots(outdir,Results);
end

Files=struct('outputDir',outdir,'summary',summaryFile,'cases',caseFile, ...
    'macro',macroFile,'P',pFile,'mat',matFile,'plots',{plotFiles});

fprintf('\nREV.4 exports written to:\n   %s\n',outdir);
end

function local_write_summary(file,Results)
fid=fopen(file,'w');
if fid<0, error('Rev4:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

C=Results.Certificate;
E=Results.Extended;

fprintf(fid,'REV.4 FINAL NUMERICAL VALIDATION SUMMARY\n');
fprintf(fid,'========================================\n\n');
fprintf(fid,'OVERALL STATUS: %s\n\n',Results.status);

fprintf(fid,'Theory synchronization\n');
fprintf(fid,'----------------------\n');
fprintf(fid,'Architecture: local flight-control macro + autonomous wind-estimator cascade\n');
fprintf(fid,'Theory file: %s\n',Results.Preflight.R.theoryFile);
fprintf(fid,'Theory SHA-256: %s\n',Results.Preflight.R.theorySHA256);
fprintf(fid,'Canonical K_omega,f: 7.0 I\n');
fprintf(fid,'certified-small maneuver scale: %.8f (validation only)\n\n',C.R.certScale);

fprintf(fid,'Local theorem-instance numerical certificate\n');
fprintf(fid,'--------------------------------------------\n');
fprintf(fid,'Sampled macro matrices: %d\n',numel(C.Family));
fprintf(fid,'Worst max Re eig(A_C,k): %.12e\n',C.maxReal);
fprintf(fid,'LMI feasp tmin: %.12e\n',C.tmin);
fprintf(fid,'Normalized common P_C min eig: %.12e\n',C.minP);
fprintf(fid,'Normalized common P_C condition number: %.12e\n',C.condP);
fprintf(fid,'Worst direct q_k: %.12e\n',C.qmin);
fprintf(fid,'Scale-invariant q_min/lambda_max(P_C): %.12e 1/s\n',C.decayProxy);
fprintf(fid,'Worst-point finite-difference q values:');
fprintf(fid,' %.12e',C.qStep);
fprintf(fid,'\n');
fprintf(fid,'Max numerical ||df_w/dx_C||_2: %.12e\n',C.maxMacroToWind);
fprintf(fid,'Analytic rho(D_M^{-1}G_M): 0 exactly\n');
fprintf(fid,'Reference-lifting min lambda/min phi: %.12f / %.12e\n',C.refMinLambda,C.refMinPhi);
fprintf(fid,'Realized local min lambda/min phi: %.12f / %.12e\n',C.realMinLambda,C.realMinPhi);
fprintf(fid,'Realized local max normalized ||X_red||: %.12e\n',C.realMaxNormX);
fprintf(fid,'Realized local max attitude tube: %.12e deg\n\n',rad2deg(C.realMaxTube));

fprintf(fid,'Extended nonlinear validation\n');
fprintf(fid,'-----------------------------\n');
for k=1:numel(E.cases)
    Q=E.cases(k);
    fprintf(fid,'%s [%s]: late Linf e_r=[%.9g %.9g %.9g] m, tube=%.9g deg, minPhi=%.9g, minLambda=%.9g\n', ...
        Q.mode,Q.role,Q.latePosLinf,Q.maxTubeDeg,Q.minPhi,Q.minLambda);
end
fprintf(fid,'\n');

fprintf(fid,'Physical residual audit (varying case)\n');
fprintf(fid,'--------------------------------------\n');
fprintf(fid,'max/late max ||d_tau^B||: %.12e / %.12e N m\n', ...
    E.varyingPhysical.maxDTau,E.varyingPhysical.lateDTau);
fprintf(fid,'max ||d_alloc||: %.12e N m\n',E.varyingPhysical.maxDAlloc);
fprintf(fid,'max omitted surface-force residual: %.12e N\n',E.varyingPhysical.maxSurfaceForce);
fprintf(fid,'max thrust lag: %.12e N\n',E.varyingPhysical.maxThrustLag);
fprintf(fid,'allocator-active samples: %d\n\n',E.varyingPhysical.allocatorActiveSamples);

fprintf(fid,'Interpretation lock\n');
fprintf(fid,'-------------------\n');
fprintf(fid,'1. The theory remains primary; the common-P LMI is a sampled T-2 instance witness only.\n');
fprintf(fid,'2. Only certified_small is used to instantiate the local macro-cascade theorem conditions.\n');
fprintf(fid,'3. varying/noisy/wrap/constant_cold are extended nonlinear robustness/performance tests.\n');
fprintf(fid,'4. Old eight-block scalar small-gain rho tests are historical diagnostics, not the final certificate.\n');
fprintf(fid,'5. The legacy 0.0405-Nm VI-D residual level is not a global full-loop physical disturbance bound.\n');
end

function local_write_cases(file,Cases)
fid=fopen(file,'w');
if fid<0, error('Rev4:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

fprintf(fid,['mode,role,late_er_N_m,late_er_E_m,late_er_D_m,' ...
    'rms_er_N_m,rms_er_E_m,rms_er_D_m,late_pos_norm_max_m,' ...
    'late_wind_error_max_mps,min_lambda,min_phi,' ...
    'thetaF_max_deg,thetaR_max_deg,tube_max_deg,' ...
    'V_min_mps,V_max_mps,gamma_abs_max_deg,alpha_min_deg,alpha_max_deg,' ...
    'beta_abs_max_deg,kappa_min,T_min_N,T_max_N,max_thrust_lag_N,max_thrust_rate_Nps,' ...
    'delta1_max_deg,delta2_max_deg,delta3_max_deg,' ...
    'delta1_rate_max_degps,delta2_rate_max_degps,delta3_rate_max_degps,' ...
    'allocator_active_samples,max_dTau_Nm,late_dTau_Nm,max_dAlloc_Nm,max_surface_force_N,' ...
    'max_slide_identity,max_wind_kinematic,max_qnorm_error,max_qrnorm_error,' ...
    'proxy_bound_N_m,proxy_bound_E_m,proxy_bound_D_m,' ...
    'proxy_ratio_N,proxy_ratio_E,proxy_ratio_D\n']);

for k=1:numel(Cases)
    Q=Cases(k);
    fprintf(fid,'%s,%s',Q.mode,Q.role);
    vals=[Q.latePosLinf Q.latePosRMS Q.latePosNormMax Q.lateWindErrorMax ...
          Q.minLambda Q.minPhi Q.maxThetaFDeg Q.maxThetaRDeg Q.maxTubeDeg ...
          Q.Vmin Q.Vmax Q.gammaAbsMaxDeg Q.alphaMinDeg Q.alphaMaxDeg ...
          Q.betaAbsMaxDeg Q.kappaMin Q.Tmin Q.Tmax Q.maxThrustLag Q.maxThrustRate ...
          Q.maxDeltaActDeg Q.maxDeltaRateDegS Q.allocatorActiveSamples ...
          Q.maxDTau Q.lateDTau Q.maxDAlloc Q.maxSurfaceForce ...
          Q.maxSlidingIdentityResidual Q.maxWindKinematicResidual ...
          Q.maxPlantQuatNormError Q.maxFilterQuatNormError ...
          Q.positionBoundProxy Q.positionProxyRatio];
    fprintf(fid,',%.16e',vals);
    fprintf(fid,'\n');
end
end

function local_write_macro(file,Fam)
fid=fopen(file,'w');
if fid<0, error('Rev4:Export','Cannot create %s',file); end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>

fprintf(fid,'index,kind,time_s,max_real_eig_Ac,q_common_P,lambda,min_phi,normX,allocator_active,tube_rad\n');
for k=1:numel(Fam)
    F=Fam(k);
    fprintf(fid,'%d,%s,%.9f,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e\n', ...
        k,F.kind,F.t,F.maxReal,F.q,F.lambda,F.minPhi,F.normX,F.allocatorActive,F.tube);
end
end

function files=local_make_plots(outdir,Results)
files={};

V=Results.Regression.CL.varying;
N=Results.Regression.CL.noisy;
C=Results.Certificate;

% 1. Position errors.
f=figure('Visible','off');
subplot(3,1,1); plot(V.t,V.Y.er(:,1),'-',N.t,N.Y.er(:,1),'--'); grid on; ylabel('e_N [m]');
subplot(3,1,2); plot(V.t,V.Y.er(:,2),'-',N.t,N.Y.er(:,2),'--'); grid on; ylabel('e_E [m]');
subplot(3,1,3); plot(V.t,V.Y.er(:,3),'-',N.t,N.Y.er(:,3),'--'); grid on; ylabel('e_D [m]'); xlabel('t [s]');
p=fullfile(outdir,'REV4_Fig1_PositionErrors.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>

% 2. Wind estimation errors.
f=figure('Visible','off');
plot(V.t,sqrt(sum(V.Y.tildew.^2,2)),'-',N.t,sqrt(sum(N.Y.tildew.^2,2)),'--');
grid on; xlabel('t [s]'); ylabel('||w - w_hat|| [m/s]','Interpreter','none');
p=fullfile(outdir,'REV4_Fig2_WindEstimation.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>

% 3. Governor feasibility.
f=figure('Visible','off');
subplot(2,1,1); plot(V.t,V.Y.lambda,'-',N.t,N.Y.lambda,'--'); grid on; ylabel('\lambda^*');
subplot(2,1,2); plot(V.t,V.Y.minPhi,'-',N.t,N.Y.minPhi,'--'); grid on; ylabel('min \phi_j'); xlabel('t [s]');
p=fullfile(outdir,'REV4_Fig3_GovernorMargins.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>

% 4. Attitude tubes.
f=figure('Visible','off');
subplot(2,1,1);
plot(V.t,rad2deg(V.Y.thetaF),'-',V.t,rad2deg(V.Y.thetaR),'--',V.t,rad2deg(V.Y.thetaF+V.Y.thetaR),':');
grid on; ylabel('varying [deg]');
subplot(2,1,2);
plot(N.t,rad2deg(N.Y.thetaF),'-',N.t,rad2deg(N.Y.thetaR),'--',N.t,rad2deg(N.Y.thetaF+N.Y.thetaR),':');
grid on; ylabel('noisy [deg]'); xlabel('t [s]');
p=fullfile(outdir,'REV4_Fig4_AttitudeTubes.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>

% 5. Physical actuator realization, varying case.
f=figure('Visible','off');
subplot(2,1,1); plot(V.t,rad2deg(V.Y.deltaAct)); grid on; ylabel('\delta_{act} [deg]');
subplot(2,1,2); plot(V.t,V.Y.Tact,'-',V.t,V.Y.Tcmd,'--'); grid on; ylabel('T [N]'); xlabel('t [s]');
p=fullfile(outdir,'REV4_Fig5_Actuators.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>

% 6. Macro certificate.
maxR=zeros(numel(C.Family),1);
q=zeros(numel(C.Family),1);
for k=1:numel(C.Family)
    maxR(k)=C.Family(k).maxReal;
    q(k)=C.Family(k).q;
end
f=figure('Visible','off');
subplot(2,1,1); plot(1:numel(maxR),maxR,'o-'); grid on; ylabel('max Re eig(A_C)');
subplot(2,1,2); plot(1:numel(q),q,'o-'); grid on; ylabel('min eig(Q_C)'); xlabel('sample index');
p=fullfile(outdir,'REV4_Fig6_MacroCertificate.png'); print(f,p,'-dpng','-r300'); close(f); files{end+1}=p; %#ok<AGROW>
end
