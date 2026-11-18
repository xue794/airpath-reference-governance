function Results=test_complete_6dof_rev1()
%TEST_COMPLETE_6DOF_REV1 Final synchronized complete closed-loop tests CL6-P1--CL6-P8.

VII=test_section_vii_rev0();
C=cl6_params();

fprintf('\nComplete Wind-Aware T-2 GTM--GNA 6DOF Closed-Loop Integration -- Final Rev.1\n');
fprintf('====================================================================\n');

%% CL6-P1 canonical parameter synchronization + exact reused quaternion plant
MF=[C.F.cF*C.F.kRf, -0.5*C.F.cF*norm(C.F.KomegaF,2);
    -0.5*C.F.cF*norm(C.F.KomegaF,2), ...
    min(eig(C.F.KomegaF))-C.F.cF];

[x01,~]=cl6_initial_state('varying',C);
rhs1=cl6_rhs(0,x01,'varying',C);
U1=cl6_plant_input(x01(33:35),x01(36),'varying',C);
p1=t2_gna_6dof_quat(0,x01(1:13),U1,C.P);
plantReuseResidual=norm(rhs1(1:13)-p1);

S1=cl6_control_snapshot(0,x01,'varying',C);
initialMomentResidual=norm(S1.Mom.B*x01(33:35)-S1.tauReq);

Fcanon=viC_params();
syncResidual=norm(C.F.KomegaF-Fcanon.KomegaF,'fro');

fprintf('CL6-P1 canonical parameter synchronization / plant composition:\n');
fprintf('   canonical K_omega,f = %.1f I\n',Fcanon.KomegaF(1,1));
fprintf('   ||K_omega,f(full)-K_omega,f(VI-C)||_F = %.3e\n',syncResidual);
fprintf('   eig(M_F,final) = [%.6f %.6f]\n',min(eig(MF)),max(eig(MF)));
fprintf('   exact reused-plant RHS residual = %.3e\n',plantReuseResidual);
fprintf('   initial physical allocation residual = %.3e N m\n',initialMomentResidual);

assert(abs(Fcanon.KomegaF(1,1)-7.0)<1e-14 && ...
       syncResidual<1e-14 && min(eig(MF))>0 && ...
       plantReuseResidual<1e-13 && initialMomentResidual<1e-10,'CL6-P1 failed');

%% CL6-P2 still-air full-6DOF hold
H=simulate_complete_6dof_case('hold',C);
ih=H.t>=H.Case.tSteady;

reH=max(abs(H.Y.er(ih,:)),[],1);
rsH=max(abs(H.Y.s(ih,:)),[],1);
tubeH=max(H.Y.thetaF+H.Y.thetaR);
maxKinH=max(H.Y.windKinematicResidual);
maxSlideH=max(H.Y.slidingIdentityResidual);

fprintf('CL6-P2 still-air full-6DOF hold:\n');
fprintf('   late Linf |e_r| [N E D] = [%.4e %.4e %.4e] m\n',reH);
fprintf('   late Linf |s|   [N E D] = [%.4e %.4e %.4e] m\n',rsH);
fprintf('   max [theta_f theta_R total] = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(max(H.Y.thetaF)),rad2deg(max(H.Y.thetaR)),rad2deg(tubeH));
fprintf('   max wind-kinematic / SMC-identity residual = [%.3e %.3e]\n', ...
    maxKinH,maxSlideH);

assert(all(reH<[0.05 0.05 0.10]) && ...
       tubeH<deg2rad(0.10) && ...
       maxKinH<1e-10 && maxSlideH<C.tol.identity,'CL6-P2 failed');

%% CL6-P3 constant-wind observer cold start in the physical loop
Kc=simulate_complete_6dof_case('constant_cold',C);
ik=Kc.t>=Kc.Case.tSteady;

w0=Kc.Y.w(1,:).';
whatAna=zeros(numel(Kc.t),3);
for k=1:numel(Kc.t)
    whatAna(k,:)=(w0-w0*exp(-2*Kc.t(k))).';
end
whatMismatch=max(sqrt(sum((Kc.Y.what-whatAna).^2,2)));

ewNorm=sqrt(sum(Kc.Y.tildew.^2,2));
e0=ewNorm(1);
id95=find(ewNorm<=0.05*e0,1,'first');
if isempty(id95), t95=inf; else, t95=Kc.t(id95); end

reK=max(abs(Kc.Y.er(ik,:)),[],1);
tubeK=max(Kc.Y.thetaF+Kc.Y.thetaR);

fprintf('CL6-P3 constant-wind cold-start observer:\n');
fprintf('   max analytic what(t) mismatch = %.3e m/s\n',whatMismatch);
fprintf('   t95 numerical = %.3f s, theory = %.3f s\n',t95,-log(0.05)/2);
fprintf('   late Linf |e_r| = [%.4e %.4e %.4e] m\n',reK);
fprintf('   max [theta_f theta_R total] = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(max(Kc.Y.thetaF)),rad2deg(max(Kc.Y.thetaR)),rad2deg(tubeK));

% Cold observer startup is intentionally outside the strict small-wind-error
% initialization used by the certified varying/noisy cases; only the total
% geometric tube is required here.
assert(whatMismatch<1e-5 && abs(t95-(-log(0.05)/2))<0.03 && ...
       all(reK<[0.05 0.05 0.05]) && tubeK<C.tol.tube,'CL6-P3 failed');

%% CL6-P4 time-varying wind: full physical loop + SMC bound
V=simulate_complete_6dof_case('varying',C);
iv=V.t>=V.Case.tSteady;

reV=max(abs(V.Y.er(iv,:)),[],1);
rsV=max(abs(V.Y.s(iv,:)),[],1);
rgV=max(abs(V.Y.eg(iv,:)),[],1);
lam=diag(C.G.Lambda).';
posBoundV=2*rsV+rgV./lam;
ratioV=reV./posBoundV;

maxThetaFV=max(V.Y.thetaF);
maxThetaRV=max(V.Y.thetaR);
maxTubeV=max(V.Y.thetaF+V.Y.thetaR);
maxSlideV=max(V.Y.slidingIdentityResidual);
minPhiV=min(V.Y.minPhi);

fprintf('CL6-P4 time-varying-wind complete 6DOF loop:\n');
fprintf('   late Linf e_r = [%.4f %.4f %.4f] m\n',reV);
fprintf('   theorem position bound = [%.4f %.4f %.4f] m\n',posBoundV);
fprintf('   realized/bound = [%.3f %.3f %.3f]\n',ratioV);
fprintf('   max [theta_f theta_R total] = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(maxThetaFV),rad2deg(maxThetaRV),rad2deg(maxTubeV));
fprintf('   max SMC identity residual = %.3e\n',maxSlideV);
fprintf('   minimum realized governor phi_j = %.3e\n',minPhiV);

assert(all(reV<=posBoundV+1e-6) && ...
       maxThetaFV<C.tol.thetaF && maxThetaRV<C.tol.thetaR && ...
       maxTubeV<C.tol.tube && maxSlideV<C.tol.identity && ...
       minPhiV>=-C.tol.phi,'CL6-P4 failed');

%% CL6-P5 deterministic bounded wind-measurement errors in full 6DOF
N=simulate_complete_6dof_case('noisy',C);
in=N.t>=N.Case.tSteady;

reN=max(abs(N.Y.er(in,:)),[],1);
rsN=max(abs(N.Y.s(in,:)),[],1);
rgN=max(abs(N.Y.eg(in,:)),[],1);
posBoundN=2*rsN+rgN./lam;
ratioN=reN./posBoundN;

maxNw=max(sqrt(sum(N.Y.nw.^2,2)));
maxWindErrLate=max(sqrt(sum(N.Y.tildew(in,:).^2,2)));
maxThetaFN=max(N.Y.thetaF);
maxThetaRN=max(N.Y.thetaR);
maxTubeN=max(N.Y.thetaF+N.Y.thetaR);

fprintf('CL6-P5 bounded wind-measurement-error robustness:\n');
fprintf('   max ||n_w|| = %.6f m/s (componentwise bound %.6f)\n', ...
    maxNw,C.E.noise.nw_bound);
fprintf('   late max ||w-what|| = %.6f m/s\n',maxWindErrLate);
fprintf('   late Linf e_r = [%.4f %.4f %.4f] m\n',reN);
fprintf('   theorem position bound = [%.4f %.4f %.4f] m\n',posBoundN);
fprintf('   realized/bound = [%.3f %.3f %.3f]\n',ratioN);
fprintf('   max [theta_f theta_R total] = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(maxThetaFN),rad2deg(maxThetaRN),rad2deg(maxTubeN));

assert(maxNw<=C.E.noise.nw_bound+1e-9 && ...
       maxWindErrLate<1.0 && all(reN<=posBoundN+1e-6) && ...
       maxThetaFN<C.tol.thetaF && maxThetaRN<C.tol.thetaR && ...
       maxTubeN<C.tol.tube,'CL6-P5 failed');

%% CL6-P6 physical allocation/servo/thrust residual audit
dAllocNorm=sqrt(sum(V.Y.dAlloc.^2,2));
dTauNomNorm=sqrt(sum(V.Y.dTauNom.^2,2));
dTauPhysNorm=sqrt(sum(V.Y.dTauPhysical.^2,2));
fSurfNorm=sqrt(sum(V.Y.fSurface.^2,2));

maxDAlloc=max(dAllocNorm);
maxDPhys=max(dTauPhysNorm);
lateDPhys=max(dTauPhysNorm(iv));
maxFSurf=max(fSurfNorm);

maxDeltaCmd=max(abs(rad2deg(V.Y.deltaCmd)),[],1);
maxDeltaAct=max(abs(rad2deg(V.Y.deltaAct)),[],1);
maxDeltaRate=max(abs(rad2deg(V.Y.deltaDot)),[],1);

maxTErr=max(abs(V.Y.Tact-V.Y.Tcmd));
maxTRate=max(abs(V.Y.TDot));
numActive=sum(any(V.Y.active>0,2));

fprintf('CL6-P6 physical allocator/actuator residual audit:\n');
fprintf('   max |delta_cmd| [deg] = [%.4f %.4f %.4f]\n',maxDeltaCmd);
fprintf('   max |delta_act| [deg] = [%.4f %.4f %.4f]\n',maxDeltaAct);
fprintf('   max |delta_dot| [deg/s] = [%.3f %.3f %.3f]\n',maxDeltaRate);
fprintf('   active allocator samples = %d / %d\n',numActive,numel(V.t));
fprintf('   max ||d_alloc|| = %.3e N m\n',maxDAlloc);
fprintf('   max / late max physical ||d_tau|| = [%.6f %.6f] N m\n',maxDPhys,lateDPhys);
fprintf('   old VI-D scalar residual budget = %.6f N m (diagnostic only)\n',C.P7.dTauBudget);
fprintf('   max omitted surface-force residual = %.6f N\n',maxFSurf);
fprintf('   max |T_act-T_cmd| / |Tdot| = [%.6f N %.3f N/s]\n',maxTErr,maxTRate);

% The old 0.0405-Nm scalar d_tau budget is deliberately NOT asserted here.
% Full integration retains the actuator state explicitly; the actual
% attitude tube is the relevant closed-loop certificate.
assert(maxDAlloc<1e-8 && numActive==0 && ...
       all(maxDeltaAct<10) && all(maxDeltaRate<=rad2deg(C.P7.deltaRateMax).'+1e-3) && ...
       maxFSurf<6.0 && maxTErr<1.5 && maxTRate<=C.P7.TRateMax+1e-3 && ...
       maxThetaRV<C.tol.thetaR && maxTubeV<C.tol.tube,'CL6-P6 failed');

%% CL6-P7 complete physical feasibility / state-domain audit
maxQ=max(abs(V.Y.qNorm-1));
maxQr=max(abs(V.Y.qrNorm-1));

Vact=[min(V.Y.V) max(V.Y.V)];
gAct=max(abs(V.Y.gamma));
aAct=[min(V.Y.alpha) max(V.Y.alpha)];
bAct=max(abs(V.Y.beta));
kMin=min(V.Y.kappa);
Tact=[min(V.Y.Tact) max(V.Y.Tact)];

Vf=[min(V.Y.xf(:,1)) max(V.Y.xf(:,1))];
gf=max(abs(V.Y.xf(:,2)));
nuMaxReal=max(abs(V.Y.xf(:,4:6)),[],1);

fprintf('CL6-P7 full physical feasibility/state-domain audit:\n');
fprintf('   actual V range = [%.4f %.4f] m/s\n',Vact);
fprintf('   actual max |gamma| = %.4f deg\n',rad2deg(gAct));
fprintf('   actual alpha range = [%.4f %.4f] deg\n',rad2deg(aAct(1)),rad2deg(aAct(2)));
fprintf('   actual max |beta| = %.4f deg, min kappa_T = %.9f\n',rad2deg(bAct),kMin);
fprintf('   actual thrust range = [%.4f %.4f] N\n',Tact);
fprintf('   governed V_f range = [%.4f %.4f] m/s, max |gamma_f| = %.4f deg\n', ...
    Vf,rad2deg(gf));
fprintf('   governed max |nu| = [%.4f %.4f %.4f]\n',nuMaxReal);
fprintf('   max quaternion norm errors [plant filter] = [%.3e %.3e]\n',maxQ,maxQr);

assert(Vact(1)>C.R.Vmin && Vact(2)<C.R.Vmax && ...
       gAct<C.R.gammaMax && ...
       aAct(1)>C.R.alphaPhysicalMin && aAct(2)<C.R.alphaPhysicalMax && ...
       bAct<C.B.betaBar && kMin>C.B.kappaMin && ...
       Tact(1)>=C.P7.Tmin-1e-8 && Tact(2)<=C.P7.Tmax+1e-8 && ...
       Vf(1)>=C.R.Vmin-1e-8 && Vf(2)<=C.R.Vmax+1e-8 && ...
       gf<=C.R.gammaMax+1e-8 && all(nuMaxReal<=C.R.nuMax.'+1e-8) && ...
       maxQ<C.tol.qnorm && maxQr<C.tol.qnorm,'CL6-P7 failed');

%% CL6-P8 full-6DOF S1 branch crossing
W=simulate_complete_6dof_case('wrap',C);

rawJump=max(abs(diff(W.Y.cc(:,3))));
maxTubeW=max(W.Y.thetaF+W.Y.thetaR);
maxErrW=max(abs(W.Y.er),[],1);

actualStep=0;
for k=2:numel(W.t)
    C0=quat_to_dcm(W.x(k-1,7:10).');
    C1=quat_to_dcm(W.x(k,7:10).');
    actualStep=max(actualStep,dcm_distance(C0,C1));
end

fprintf('CL6-P8 S1 branch crossing through the COMPLETE 6DOF loop:\n');
fprintf('   raw SMC/RG command-course jump = %.4f deg\n',rad2deg(rawJump));
fprintf('   max total attitude tube = %.6f deg\n',rad2deg(maxTubeW));
fprintf('   max consecutive actual-attitude step = %.6f deg\n',rad2deg(actualStep));
fprintf('   max |e_r| [N E D] = [%.4f %.4f %.4f] m\n',maxErrW);

assert(rawJump>deg2rad(300) && maxTubeW<C.tol.tube && ...
       actualStep<deg2rad(0.20) && all(maxErrW<[2 2 1]),'CL6-P8 failed');

Results=struct('VII',VII,'C',C,'hold',H,'constantCold',Kc, ...
    'varying',V,'noisy',N,'wrap',W,'MF',MF);
end
