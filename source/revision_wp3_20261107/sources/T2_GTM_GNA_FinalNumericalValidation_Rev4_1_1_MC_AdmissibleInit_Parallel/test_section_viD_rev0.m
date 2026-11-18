function Results=test_section_viD_rev0()
%TEST_SECTION_VID_REV0 Section VI-D tests VID-P1--VID-P8.

VIC=test_section_viC_rev0();
R=VIC.R;
A=VIC.A;
B=VIC.B;
F=VIC.F;
Dp=viD_params();
Ew=wind_estimator_params();

fprintf('\nSection VI-D PHS-Based Rotational Attitude/Rate Controller Rev.0\n');
fprintf('=================================================================\n');

%% VID-P1 exact T-2 GNA rate-moment decomposition and D_omega,a positivity
Cert=viD_damping_certificate(R,Dp);

rng(51,'twister');
maxGNA=0; maxDecomp=0; maxSkew=0; maxSym=0;
for k=1:500
    V=R.Vmin+(R.Vmax-R.Vmin)*rand;
    alpha=R.alphaPhysicalMin+(R.alphaPhysicalMax-R.alphaPhysicalMin)*rand;
    beta=deg2rad(-1+2*rand);
    omega=-0.3+0.6*rand(3,1);

    va=V*[cos(alpha)*cos(beta);sin(beta);sin(alpha)*cos(beta)];
    [~,tau0]=t2_gna_aero(va,zeros(3,1),zeros(3,1),eye(3),zeros(3,1),R.h_ref,R.P);
    [~,tau1]=t2_gna_aero(va,omega,zeros(3,1),eye(3),zeros(3,1),R.h_ref,R.P);

    M=viD_rate_moment_decomp(V,alpha,R);

    maxGNA=max(maxGNA,norm((tau1-tau0)-M.K*omega));
    maxDecomp=max(maxDecomp,norm(M.K+M.D-M.S,'fro'));
    maxSkew=max(maxSkew,norm(M.S+M.S.','fro'));
    maxSym=max(maxSym,norm(M.D-M.D.','fro'));
end

fprintf('VID-P1 T-2 GNA rotational-rate PHS decomposition:\n');
fprintf('   max GNA rate-moment reconstruction residual = %.3e N m\n',maxGNA);
fprintf('   max K+D-S residual = %.3e\n',maxDecomp);
fprintf('   max skew/symmetric structural residuals = [%.3e %.3e]\n',maxSkew,maxSym);
fprintf('   grid min eig(D_omega,a) = %.6f N m s/rad\n',Cert.minEigD);
fprintf('   grid max ||D_omega,a||_2 = %.6f N m s/rad\n',Cert.maxNormD);
assert(maxGNA<1e-12 && maxDecomp<1e-12 && maxSkew<1e-12 && ...
       maxSym<1e-12 && Cert.minEigD>0,'VID-P1 failed');

%% VID-P2 rotational error kinematics and exact virtual-torque reduction
rng(52,'twister');
maxQdot=0; maxEdot=0; maxPsiDot=0; maxDyn=0;
for k=1:500
    Cr=rotvec_to_dcm(randn(3,1)*0.4);
    ax=randn(3,1); ax=ax/norm(ax);
    QR=rotvec_to_dcm(ax*deg2rad(5*rand));
    C=Cr*QR;

    omegaR=-0.10+0.20*rand(3,1);
    omegaRDot=-0.30+0.60*rand(3,1);

    Err0=viD_error(Cr,omegaR,C,zeros(3,1));
    omega=Err0.A*omegaR+(-0.05+0.10*rand(3,1));
    Err=viD_error(Cr,omegaR,C,omega);

    V=32+10*rand;
    alpha=deg2rad(1+7*rand);
    Mom=viD_rate_moment_decomp(V,alpha,R);
    dTau=-0.05+0.10*rand(3,1);

    Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot,C,omega,Mom.D,R,Dp);
    omegaDot=R.P.I\(-cross(omega,R.P.I*omega)-Mom.D*omega+Ctl.tauC+dTau);

    Qdot=-viC_hat(omegaR)*Err.Q+Err.Q*viC_hat(omega);
    QdotTheory=Err.Q*viC_hat(Err.eOmega);
    eDot=0.5*viC_vee(Qdot-Qdot.');
    eDotTheory=viC_E(Err.Q)*Err.eOmega;
    psiDot=-0.5*trace(Qdot);

    eOmegaDot=omegaDot + viC_hat(omega)*(Err.A*omegaR) - Err.A*omegaRDot;
    eOmegaDotTheory=R.P.I\(-Dp.kR*Err.eR-(Dp.KomegaC+Mom.D)*Err.eOmega+dTau);

    maxQdot=max(maxQdot,norm(Qdot-QdotTheory,'fro'));
    maxEdot=max(maxEdot,norm(eDot-eDotTheory));
    maxPsiDot=max(maxPsiDot,abs(psiDot-Err.eR.'*Err.eOmega));
    maxDyn=max(maxDyn,norm(eOmegaDot-eOmegaDotTheory));
end

fprintf('VID-P2 rotational error / virtual-torque identities:\n');
fprintf('   max Q_R-dot residual       = %.3e\n',maxQdot);
fprintf('   max e_R-dot residual       = %.3e\n',maxEdot);
fprintf('   max Psi_R-dot residual     = %.3e\n',maxPsiDot);
fprintf('   max e_omega dynamics residual = %.3e rad/s^2\n',maxDyn);
assert(max([maxQdot maxEdot maxPsiDot maxDyn])<1e-12,'VID-P2 failed');

%% VID-P3 strict local rotational storage / ISS inequality
I=R.P.I;
Imin=min(eig(I)); Imax=norm(I,2);

KdMin=Dp.KomegaC(1,1)+Cert.minEigD;
KdMax=norm(Dp.KomegaC,2)+Cert.maxNormD;

PR=[0.5*Dp.kR, -0.5*Dp.cR*Imax;
    -0.5*Dp.cR*Imax, 0.5*Imin];

MR=[Dp.cR*Dp.kR, -0.5*Dp.cR*KdMax;
    -0.5*Dp.cR*KdMax, KdMin-Dp.cR*Imax];

rng(53,'twister');
minStorageGap=inf; maxDerivViolation=-inf;
for k=1:5000
    ax=randn(3,1); ax=ax/norm(ax);
    Q=rotvec_to_dcm(ax*(Dp.thetaRMax*rand));
    eR=0.5*viC_vee(Q-Q.');
    Psi=0.5*trace(eye(3)-Q);
    eOmega=-0.05+0.10*rand(3,1);

    V=R.Vmin+(R.Vmax-R.Vmin)*rand;
    alpha=R.alphaPhysicalMin+(R.alphaPhysicalMax-R.alphaPhysicalMin)*rand;
    Mom=viD_rate_moment_decomp(V,alpha,R);
    Kd=Dp.KomegaC+Mom.D;

    dTau=0.04*(2*rand(3,1)-1);
    eDot=viC_E(Q)*eOmega;
    eOmegaDot=I\(-Dp.kR*eR-Kd*eOmega+dTau);

    VR=0.5*eOmega.'*I*eOmega + Dp.kR*Psi + Dp.cR*(eR.'*I*eOmega);
    Vdot=eOmega.'*I*eOmegaDot + Dp.kR*(eR.'*eOmega) + ...
         Dp.cR*(eDot.'*I*eOmega + eR.'*I*eOmegaDot);

    qn=[norm(eR);norm(eOmega)];
    lower=qn.'*PR*qn;
    upper=-qn.'*MR*qn + (qn(2)+Dp.cR*qn(1))*norm(dTau);

    minStorageGap=min(minStorageGap,VR-lower);
    maxDerivViolation=max(maxDerivViolation,Vdot-upper);
end

fprintf('VID-P3 rotational storage / ISS certificate:\n');
fprintf('   eig(P_R lower matrix) = [%.6f %.6f]\n',min(eig(PR)),max(eig(PR)));
fprintf('   eig(M_R)              = [%.6f %.6f]\n',min(eig(MR)),max(eig(MR)));
fprintf('   sampled min storage lower-bound gap = %.3e\n',minStorageGap);
fprintf('   max derivative-bound violation       = %.3e\n',maxDerivViolation);
assert(min(eig(PR))>0 && min(eig(MR))>0 && ...
       minStorageGap>-1e-12 && maxDerivViolation<1e-10,'VID-P3 failed');

%% VID-P4 fixed-reference convergence inside theta_R tube
Cr=viB_command_attitude(deg2rad(2),deg2rad(15),deg2rad(5),deg2rad(5));
QR0=rotvec_to_dcm(Dp.initialAxis*Dp.constTheta0);
C0=Cr*QR0;
q0=viC_dcm_to_quat(C0);
omega0=Dp.constEOmega0;

tg=linspace(0,Dp.Tconst,1501);
opts=odeset('RelTol',Dp.RelTol,'AbsTol',Dp.AbsTol,'MaxStep',0.002);
[t4,x4]=ode45(@(tt,xx)viD_constant_rhs(tt,xx,Cr,R.Vb,deg2rad(5),R,F,Dp), ...
              tg,[q0;omega0],opts);

theta4=zeros(numel(t4),1);
for k=1:numel(t4)
    theta4(k)=dcm_distance(Cr,quat_to_dcm(x4(k,1:4).'));
end

maxTheta4=max(theta4);
finalTheta4=theta4(end);
finalOmega4=norm(x4(end,5:7));
maxQnorm4=max(abs(sqrt(sum(x4(:,1:4).^2,2))-1));

fprintf('VID-P4 constant-reference rotational convergence:\n');
fprintf('   max theta_R = %.6f deg (tube %.6f deg)\n',rad2deg(maxTheta4),rad2deg(Dp.thetaRMax));
fprintf('   final theta_R = %.3e rad\n',finalTheta4);
fprintf('   final |omega| = %.3e rad/s\n',finalOmega4);
fprintf('   max quaternion-norm error = %.3e\n',maxQnorm4);
assert(maxTheta4<Dp.thetaRMax && finalTheta4<5e-8 && ...
       finalOmega4<1e-8 && maxQnorm4<1e-8,'VID-P4 failed');

%% VID-P5 aggressive moving reference, ideal torque port, no residual
S5=simulate_viD_case('aggressive',Dp.Taggressive,R,A,F,Dp,Ew);
maxThetaF5=max(S5.thetaF);
maxThetaR5=max(S5.thetaR);
maxTube5=max(S5.thetaF+S5.thetaR);
lateR5=max(S5.thetaR(S5.t>=Dp.tSteady));
maxEOmega5=max(sqrt(sum(S5.eOmega.^2,2)));
maxTauAxis5=max(abs(S5.tauC),[],1);

fprintf('VID-P5 aggressive moving-reference rotational tracking:\n');
fprintf('   max [theta_f theta_R theta_f+theta_R] = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(maxThetaF5),rad2deg(maxThetaR5),rad2deg(maxTube5));
fprintf('   late max theta_R = %.6e deg\n',rad2deg(lateR5));
fprintf('   max |e_omega| = %.6f rad/s\n',maxEOmega5);
fprintf('   max |tau_c| per axis = [%.4f %.4f %.4f] N m\n',maxTauAxis5);
assert(maxThetaF5<F.thetaFMax && maxThetaR5<Dp.thetaRMax && ...
       maxTube5<B.thetaTube && lateR5<1e-6,'VID-P5 failed');

%% VID-P6 bounded torque residual: ISS/UUB and tube retention
S6=simulate_viD_case('aggressive_disturbed',Dp.Taggressive,R,A,F,Dp,Ew);
idx6=S6.t>=Dp.tSteady;
maxThetaR6=max(S6.thetaR);
lateR6=max(S6.thetaR(idx6));
rmsR6=sqrt(mean(S6.thetaR(idx6).^2));
maxTube6=max(S6.thetaF+S6.thetaR);
dNorm=sqrt(sum(S6.dTau.^2,2));
dBar=max(dNorm);

lambdaMR=min(eig(MR));
qUltimate=sqrt(1+Dp.cR^2)*dBar/lambdaMR;
thetaUltimate=asin(min(1,qUltimate));

fprintf('VID-P6 bounded rotational residual ISS/UUB:\n');
fprintf('   max ||d_tau|| = %.6f N m\n',dBar);
fprintf('   realized max / late max / late RMS theta_R = [%.6f %.6f %.6f] deg\n', ...
    rad2deg(maxThetaR6),rad2deg(lateR6),rad2deg(rmsR6));
fprintf('   conservative storage q-radius = %.6e, theta equivalent %.6f deg\n', ...
    qUltimate,rad2deg(thetaUltimate));
fprintf('   max theta_f+theta_R = %.6f deg\n',rad2deg(maxTube6));
assert(maxThetaR6<Dp.thetaRMax && maxTube6<B.thetaTube && ...
       thetaUltimate<Dp.thetaRMax,'VID-P6 failed');

%% VID-P7 actual alpha/beta/kappa/thrust transfer with dynamic Q_f and Q_R
amin=min(S6.alpha); amax=max(S6.alpha);
bmax=max(abs(S6.beta));
kmin=min(S6.kappa);
Tmin=min(S6.Tphys); Tmax=max(S6.Tphys);
thetaCmax=max(S6.thetaC);
triGap=max(S6.thetaC-(S6.thetaF+S6.thetaR));

fprintf('VID-P7 full dynamic attitude-geometry feasibility transfer:\n');
fprintf('   max theta_c = %.6f deg, max triangle residual = %.3e rad\n', ...
    rad2deg(thetaCmax),triGap);
fprintf('   actual alpha range = [%.6f %.6f] deg\n',rad2deg(amin),rad2deg(amax));
fprintf('   max |beta| = %.6f deg, min kappa_T = %.9f\n',rad2deg(bmax),kmin);
fprintf('   physical thrust range = [%.6f %.6f] N\n',Tmin,Tmax);
assert(thetaCmax<B.thetaTube && triGap<1e-8 && ...
       amin>R.alphaPhysicalMin && amax<R.alphaPhysicalMax && ...
       bmax<B.betaBar && kmin>B.kappaMin && ...
       Tmin>R.Tmin+R.DeltaT && Tmax<R.Tmax-R.DeltaT,'VID-P7 failed');

%% VID-P8 topology / wind robustness through actual rotational loop
S8w=simulate_viD_case('wrap',Dp.Twrap,R,A,F,Dp,Ew);
S8v=simulate_viD_case('varying',Dp.Tvarying,R,A,F,Dp,Ew);

rawChi=zeros(numel(S8w.t),1);
for k=1:numel(S8w.t)
    [vc,~]=rg_raw_test_command(S8w.t(k),'wrap',R);
    [cc,~]=rg_raw_command_chart(vc,R);
    rawChi(k)=cc(3);
end
rawJump=max(abs(diff(rawChi)));

Cstep=0;
for k=2:numel(S8w.t)
    C0=quat_to_dcm(S8w.x(k-1,8:11).');
    C1=quat_to_dcm(S8w.x(k,8:11).');
    Cstep=max(Cstep,dcm_distance(C0,C1));
end

maxWrapR=max(S8w.thetaR);
maxWrapTube=max(S8w.thetaF+S8w.thetaR);
maxVarR=max(S8v.thetaR);
maxVarTube=max(S8v.thetaF+S8v.thetaR);

fprintf('VID-P8 S1 / time-varying-wind rotational robustness:\n');
fprintf('   raw course branch jump = %.4f deg\n',rad2deg(rawJump));
fprintf('   wrap: max theta_R / total tube = [%.6f %.6f] deg\n', ...
    rad2deg(maxWrapR),rad2deg(maxWrapTube));
fprintf('   wrap: max consecutive actual-attitude step = %.6f deg\n',rad2deg(Cstep));
fprintf('   varying wind: max theta_R / total tube = [%.6f %.6f] deg\n', ...
    rad2deg(maxVarR),rad2deg(maxVarTube));
assert(rawJump>deg2rad(300) && Cstep<deg2rad(0.2) && ...
       maxWrapR<Dp.thetaRMax && maxWrapTube<B.thetaTube && ...
       maxVarR<Dp.thetaRMax && maxVarTube<B.thetaTube,'VID-P8 failed');

Results=struct('VIC',VIC,'R',R,'A',A,'B',B,'F',F,'Dp',Dp, ...
    'VID1',struct('certificate',Cert,'maxGNA',maxGNA,'maxDecomp',maxDecomp), ...
    'VID2',struct('Qdot',maxQdot,'eDot',maxEdot,'PsiDot',maxPsiDot,'dyn',maxDyn), ...
    'VID3',struct('PR',PR,'MR',MR,'minStorageGap',minStorageGap, ...
                  'maxDerivViolation',maxDerivViolation), ...
    'VID4',struct('t',t4,'x',x4,'theta',theta4), ...
    'VID5',S5,'VID6',S6,'VID7',struct('amin',amin,'amax',amax,'bmax',bmax, ...
                  'kmin',kmin,'Tmin',Tmin,'Tmax',Tmax), ...
    'VID8',struct('wrap',S8w,'varying',S8v));
end
