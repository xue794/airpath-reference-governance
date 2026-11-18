function Results=test_section_viC_rev0()
%TEST_SECTION_VIC_REV0 Section VI-C SO(3) command-filter tests VIC-P1--P8.

VIB=test_section_viB_rev0();
R=VIB.R;
A=VIB.A;
B=VIB.B;
F=viC_params();
Ew=wind_estimator_params();

fprintf('\nSection VI-C SO(3) Command Filter -- Final Synchronized Rev.1\n');
fprintf('=======================================\n');

%% VIC-P1 analytic command derivative / no-hidden-third-derivative closure
Up=simulate_viA_case('aggressive',A.Taggressive,R,A,Ew);
idx=round(linspace(2,numel(Up.t)-1,F.Nkin));
maxAlphaDot=0; maxMuDot=0; maxOmega=0; maxOmegaMag=0;

for jj=1:numel(idx)
    k=idx(jj);
    t=Up.t(k);
    z=Up.z(k,:).';

    Kin=viC_command_kinematics(t,z,'aggressive',R,A,Ew);
    z1=z+F.fdDt*Kin.upDot;
    O1=viC_command_from_upstream(z1,R,A);

    alphaFD=(O1.alpha-Kin.C.alpha)/F.fdDt;
    muFD=rg_angle_error(O1.mu,Kin.C.mu)/F.fdDt;

    Qstep=Kin.Cc.'*O1.Cc;
    omegaFD=viC_vee(0.5*(Qstep-Qstep.'))/F.fdDt;

    maxAlphaDot=max(maxAlphaDot,abs(alphaFD-Kin.alphaDot));
    maxMuDot=max(maxMuDot,abs(muFD-Kin.muDot));
    maxOmega=max(maxOmega,norm(omegaFD-Kin.omegaC));
    maxOmegaMag=max(maxOmegaMag,norm(Kin.omegaC));
end

fprintf('VIC-P1 command-attitude derivative closure:\n');
fprintf('   max alpha_dot finite-difference error = %.3e rad/s\n',maxAlphaDot);
fprintf('   max mu_dot finite-difference error    = %.3e rad/s\n',maxMuDot);
fprintf('   max omega_c finite-difference error   = %.3e rad/s\n',maxOmega);
fprintf('   max |omega_c| on aggressive trajectory = %.6f rad/s\n',maxOmegaMag);
assert(maxAlphaDot<2e-8 && maxMuDot<3e-8 && maxOmega<3e-8,'VIC-P1 failed');

%% VIC-P2 exact filter-error kinematics
rng(41,'twister');
maxQ=0; maxPsi=0; maxE=0;
for k=1:1000
    ax=randn(3,1); ax=ax/norm(ax);
    Q=rotvec_to_dcm(ax*(deg2rad(20)*rand));
    xi=-0.2+0.4*rand(3,1);

    Qdot=Q*viC_hat(xi);
    Q1=Q+F.fdDt*Qdot;

    e0=0.5*viC_vee(Q-Q.');
    e1=0.5*viC_vee(Q1-Q1.');
    psi0=0.5*trace(eye(3)-Q);
    psi1=0.5*trace(eye(3)-Q1);

    maxQ=max(maxQ,norm(Qdot-Q*viC_hat(xi)));
    maxPsi=max(maxPsi,abs((psi1-psi0)/F.fdDt-e0.'*xi));
    maxE=max(maxE,norm((e1-e0)/F.fdDt-viC_E(Q)*xi));
end
fprintf('VIC-P2 SO(3) filter-error kinematics:\n');
fprintf('   max Qdot identity residual   = %.3e\n',maxQ);
fprintf('   max Psi_dot identity residual = %.3e\n',maxPsi);
fprintf('   max e_dot identity residual   = %.3e\n',maxE);
assert(maxQ<1e-12 && maxPsi<2e-8 && maxE<2e-8,'VIC-P2 failed');

%% VIC-P3 Lyapunov/storage sufficient-condition audit
KF=F.KomegaF;
MF=[F.cF*F.kRf, -0.5*F.cF*norm(KF,2);
    -0.5*F.cF*norm(KF,2), min(eig(KF))-F.cF];

evMF=eig(MF);
rng(42,'twister');
maxViolation=-inf; minStorage=inf;
for k=1:F.Nlyap
    ax=randn(3,1); ax=ax/norm(ax);
    Q=rotvec_to_dcm(ax*(deg2rad(30)*rand));
    omegaR=-0.2+0.4*rand(3,1);
    omegaC=-0.08+0.16*rand(3,1);

    e=0.5*viC_vee(Q-Q.');
    Psi=0.5*trace(eye(3)-Q);
    xi=omegaR-Q.'*omegaC;
    eDot=viC_E(Q)*xi;
    PsiDot=e.'*xi;
    omegaRDot=-F.kRf*e-KF*omegaR;

    VF=0.5*(omegaR.'*omegaR)+F.kRf*Psi+F.cF*(e.'*omegaR);
    Vdot=omegaR.'*omegaRDot+F.kRf*PsiDot+ ...
         F.cF*(eDot.'*omegaR+e.'*omegaRDot);

    q=[norm(e);norm(omegaR)];
    rhs=-q.'*MF*q+(F.kRf*q(1)+F.cF*q(2))*norm(omegaC);

    maxViolation=max(maxViolation,Vdot-rhs);
    minStorage=min(minStorage,VF);
end

fprintf('VIC-P3 storage / derivative inequality audit:\n');
fprintf('   c_F^2 < k_Rf: %.6f < %.6f\n',F.cF^2,F.kRf);
fprintf('   eig(M_F) = [%.6f %.6f]\n',min(evMF),max(evMF));
fprintf('   sampled minimum V_F = %.6e\n',minStorage);
fprintf('   max sampled derivative-bound violation = %.3e\n',maxViolation);
assert(F.cF^2<F.kRf && min(evMF)>0 && minStorage>0 && maxViolation<=1e-10,'VIC-P3 failed');

%% VIC-P4 constant-command asymptotic convergence
Cc=viB_command_attitude(deg2rad(2),deg2rad(15),deg2rad(5),deg2rad(5));
axis0=[1;-0.5;0.3]; axis0=axis0/norm(axis0);
Cr0=Cc*rotvec_to_dcm(axis0*F.constOffset);
qr0=viC_dcm_to_quat(Cr0);
x0=[qr0;F.constOmega0];

tg=linspace(0,F.constT,1001);
opts=odeset('RelTol',F.RelTol,'AbsTol',F.AbsTol,'MaxStep',0.002);
[t4,x4]=ode45(@(tt,xx)viC_constant_rhs(tt,xx,Cc,F),tg,x0,opts);

theta4=zeros(numel(t4),1);
for k=1:numel(t4)
    theta4(k)=dcm_distance(Cc,quat_to_dcm(x4(k,1:4).'));
end
finalTheta=theta4(end);
finalOmega=norm(x4(end,5:7));
maxTheta=max(theta4);
maxQnorm=max(abs(sqrt(sum(x4(:,1:4).^2,2))-1));

fprintf('VIC-P4 constant-command convergence:\n');
fprintf('   max theta_f = %.6f deg\n',rad2deg(maxTheta));
fprintf('   final theta_f = %.3e rad\n',finalTheta);
fprintf('   final |omega_r| = %.3e rad/s\n',finalOmega);
fprintf('   max quaternion-norm error = %.3e\n',maxQnorm);
assert(maxTheta<F.thetaFMax && finalTheta<5e-8 && finalOmega<1e-8 && maxQnorm<1e-8,'VIC-P4 failed');

%% VIC-P5 persistent moving command: bounded nonzero lag
x0=[1;0;0;0;0;0;0];
tg=linspace(0,F.spinT,4001);
opts=odeset('RelTol',F.RelTol,'AbsTol',F.AbsTol,'MaxStep',0.01);
[t5,x5]=ode45(@(tt,xx)viC_spin_rhs(tt,xx,F),tg,x0,opts);

theta5=zeros(numel(t5),1);
for k=1:numel(t5)
    Cc5=viB_Rz(F.spinRate*t5(k));
    theta5(k)=dcm_distance(Cc5,quat_to_dcm(x5(k,1:4).'));
end
iss=t5>=F.spinSteadyT;
steadyTheta=mean(theta5(iss));
steadyOmega=mean(x5(iss,7));
predTheta=asin(F.KomegaF(3,3)*F.spinRate/F.kRf);

fprintf('VIC-P5 persistent rotating-command UUB/nonzero lag:\n');
fprintf('   steady theta_f = %.6f deg\n',rad2deg(steadyTheta));
fprintf('   exact equilibrium prediction = %.6f deg\n',rad2deg(predTheta));
fprintf('   steady omega_r,z = %.9f vs omega_c = %.9f rad/s\n',steadyOmega,F.spinRate);
assert(abs(steadyTheta-predTheta)<2e-7 && ...
       abs(steadyOmega-F.spinRate)<1e-8 && steadyTheta>deg2rad(0.01),'VIC-P5 failed');

%% VIC-P6 actual aggressive Section-V/VI-A command: dynamic realization of VI-B tube
S6=simulate_viC_case('aggressive',A.Taggressive,R,A,F,Ew);
maxTheta6=max(S6.theta);
lateTheta6=max(S6.theta(S6.t>=60));
rmsLate6=sqrt(mean(S6.theta(S6.t>=60).^2));
maxOmegaC6=max(sqrt(sum(S6.omegaC.^2,2)));
maxOmegaR6=max(sqrt(sum(S6.omegaR.^2,2)));
maxXi6=max(sqrt(sum(S6.xi.^2,2)));
maxPsi6=max(S6.Psi);
maxQnorm6=max(abs(S6.qNorm-1));

fprintf('VIC-P6 aggressive moving command-filter realization:\n');
fprintf('   max theta_f = %.6f deg (tube %.6f deg)\n',rad2deg(maxTheta6),rad2deg(F.thetaFMax));
fprintf('   late max / RMS theta_f = [%.6f %.6f] deg\n',rad2deg(lateTheta6),rad2deg(rmsLate6));
fprintf('   max [|omega_c| |omega_r| |xi_f|] = [%.6f %.6f %.6f] rad/s\n', ...
    maxOmegaC6,maxOmegaR6,maxXi6);
fprintf('   max Psi_f = %.6e, max quaternion-norm error = %.3e\n',maxPsi6,maxQnorm6);
assert(maxTheta6<F.thetaFMax && maxPsi6<F.PsiFMax && maxQnorm6<1e-8,'VIC-P6 failed');

%% VIC-P7 S1 branch crossing with dynamic filter
S7=simulate_viC_case('wrap',A.Twrap,R,A,F,Ew);
rawChi=zeros(numel(S7.t),1);
for k=1:numel(S7.t)
    [vc,~]=rg_raw_test_command(S7.t(k),'wrap',R);
    [cc,~]=rg_raw_command_chart(vc,R);
    rawChi(k)=cc(3);
end
rawJump=max(abs(diff(rawChi)));
maxTheta7=max(S7.theta);
maxOmegaC7=max(sqrt(sum(S7.omegaC.^2,2)));

CrStep=0;
for k=2:numel(S7.t)
    C0=quat_to_dcm(S7.x(k-1,1:4).');
    C1=quat_to_dcm(S7.x(k,1:4).');
    CrStep=max(CrStep,dcm_distance(C0,C1));
end

fprintf('VIC-P7 S1 branch crossing through dynamic command filter:\n');
fprintf('   raw command-course jump = %.4f deg\n',rad2deg(rawJump));
fprintf('   max theta_f = %.6f deg\n',rad2deg(maxTheta7));
fprintf('   max |omega_c| = %.6f rad/s\n',maxOmegaC7);
fprintf('   max consecutive C_r geodesic step = %.6f deg\n',rad2deg(CrStep));
assert(rawJump>deg2rad(300) && maxTheta7<F.thetaFMax && CrStep<deg2rad(0.2),'VIC-P7 failed');

%% VIC-P8 time-varying-wind upstream command
S8=simulate_viC_case('varying',A.Tvarying,R,A,F,Ew);
idx8=S8.t>=A.tUUB0;
maxTheta8=max(S8.theta);
lateTheta8=max(S8.theta(idx8));
rmsTheta8=sqrt(mean(S8.theta(idx8).^2));
maxOmegaC8=max(sqrt(sum(S8.omegaC.^2,2)));

fprintf('VIC-P8 time-varying-wind command-filter UUB:\n');
fprintf('   max theta_f = %.6f deg\n',rad2deg(maxTheta8));
fprintf('   late max / RMS theta_f = [%.6f %.6f] deg\n',rad2deg(lateTheta8),rad2deg(rmsTheta8));
fprintf('   max |omega_c| = %.6f rad/s\n',maxOmegaC8);
assert(maxTheta8<F.thetaFMax && lateTheta8>1e-6,'VIC-P8 failed');

Results=struct('VIB',VIB,'R',R,'A',A,'B',B,'F',F, ...
    'VIC1',struct('alphaDotErr',maxAlphaDot,'muDotErr',maxMuDot, ...
                  'omegaCErr',maxOmega,'maxOmegaC',maxOmegaMag), ...
    'VIC2',struct('QdotErr',maxQ,'PsiDotErr',maxPsi,'eDotErr',maxE), ...
    'VIC3',struct('MF',MF,'eigMF',evMF,'minStorage',minStorage,'maxViolation',maxViolation), ...
    'VIC4',struct('t',t4,'x',x4,'theta',theta4), ...
    'VIC5',struct('t',t5,'x',x5,'theta',theta5,'prediction',predTheta), ...
    'VIC6',S6,'VIC7',S7,'VIC8',S8);
end
