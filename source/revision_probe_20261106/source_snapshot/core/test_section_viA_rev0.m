function Results = test_section_viA_rev0()
%TEST_SECTION_VIA_REV0 Section VI-A module tests VIA-P1--VIA-P8.

RG=test_reference_governor_rev0();
R=RG.R;
A=viA_params();
E=wind_estimator_params();

fprintf('\nSection VI-A Air-Path Force Inversion and Ideal Realization Rev.0\n');
fprintf('=================================================================\n');

%% VIA-P1 air-path basis and ideal force decomposition
rng(21,'twister');
maxOrth=0; maxForce=0;
for k=1:500
    gamma=deg2rad(-10+20*rand);
    chi=-pi+2*pi*rand;
    mu=deg2rad(-35+70*rand);

    B=viA_airpath_basis(gamma,chi);
    Q=[B.t B.nGamma B.nChi];
    maxOrth=max(maxOrth,norm(Q.'*Q-eye(3),'fro'));

    V=32+12*rand;
    alpha=deg2rad(2+6*rand);
    X=-20+40*rand;
    M=rg_N_map(alpha,V,X,R);
    Aa=rg_nominal_aero(V,alpha,R);
    Fp=Aa.D+X;
    ell=cos(mu)*B.nGamma+sin(mu)*B.nChi;
    f=(Fp-Aa.D)*B.t+M.value*ell;

    comp=[B.t.'*f;B.nGamma.'*f;B.nChi.'*f];
    target=[X;M.value*cos(mu);M.value*sin(mu)];
    maxForce=max(maxForce,norm(comp-target));
end
fprintf('VIA-P1 air-path basis / force decomposition:\n');
fprintf('   max orthonormality residual = %.3e\n',maxOrth);
fprintf('   max force-component residual = %.3e N\n',maxForce);
assert(maxOrth<1e-12 && maxForce<1e-10,'VIA-P1 failed');

%% VIA-P2 zero-error consistency with Section V feedforward inversion
rng(22,'twister');
maxA=0; maxT=0; maxMu=0; maxX=0; maxN=0;
for k=1:500
    xf=[R.Vb+(-0.5+rand);
        deg2rad(-1+2*rand);
        -pi+2*pi*rand;
        -0.1+0.2*rand;
        deg2rad(-0.5+rand);
        deg2rad(-1+2*rand)];

    F=rg_feedforward(xf,R);
    if ~F.rootInfo.bracket
        error('VIA-P2 sampled feedforward root lost.');
    end
    C=viA_feedback_command(xf(1:3),xf,R,A);

    maxA=max(maxA,abs(C.alpha-F.alpha));
    maxT=max(maxT,abs(C.Tstar-F.T));
    maxMu=max(maxMu,abs(rg_angle_error(C.mu,F.mu)));
    maxX=max(maxX,abs(C.X-F.X));
    maxN=max(maxN,abs(C.N-F.N));
end
fprintf('VIA-P2 zero-error feedforward consistency:\n');
fprintf('   max |alpha_c-alpha_ff| = %.3e rad\n',maxA);
fprintf('   max |T*-T_ff|          = %.3e N\n',maxT);
fprintf('   max |mu_c-mu_ff|       = %.3e rad\n',maxMu);
fprintf('   max |X_c-X_f|, |N_c-N_f| = [%.3e %.3e]\n',maxX,maxN);
assert(max([maxA maxT maxMu maxX maxN])<1e-11,'VIA-P2 failed');

%% VIA-P3 exact ideal error-dynamics identity
rng(23,'twister');
maxId=0;
for k=1:500
    xf=[R.Vb+(-0.4+0.8*rand);
        deg2rad(-1+2*rand);
        -2.5+5*rand;
        -0.1+0.2*rand;
        deg2rad(-0.6+1.2*rand);
        deg2rad(-1.2+2.4*rand)];

    y=[xf(1)+(-0.4+0.8*rand);
       xf(2)+deg2rad(-0.6+1.2*rand);
       xf(3)+deg2rad(-1+2*rand)];

    wd=-0.7+1.4*rand(3,1);
    [ydot,~]=viA_ideal_path_rhs(y,xf,wd,R,A);
    direct=ydot-xf(4:6);
    theory=viA_error_rhs_theory(y,xf,wd,R,A);
    maxId=max(maxId,norm(direct-theory));
end
fprintf('VIA-P3 exact ideal path-error dynamics:\n');
fprintf('   max identity residual = %.3e\n',maxId);
assert(maxId<1e-12,'VIA-P3 failed');

%% VIA-P4 still-air exponential realization of a moving governed reference
S4=simulate_viA_case('still',A.Tstill,R,A,E);
e0=A.e0;
ea=[e0(1)*exp(-A.kV*S4.t), ...
    e0(2)*exp(-A.kGamma*S4.t), ...
    e0(3)*exp(-A.kChi*S4.t)];
maxExp=max(max(abs(S4.e-ea)));
late=max(abs(S4.e(S4.t>=50,:)),[],1);
fprintf('VIA-P4 still-air moving-reference realization:\n');
fprintf('   max analytic error-trajectory mismatch = %.3e\n',maxExp);
fprintf('   late Linf [eV eGamma eChi] = [%.3e %.3e %.3e]\n',late);
assert(maxExp<5e-7 && max(late)<1e-8,'VIA-P4 failed');

%% VIA-P5 constant-wind magnitude invariance
S5=simulate_viA_case('constant',A.Tconstant,R,A,E);
maxState=max(max(abs(S4.z-S5.z)));
maxCmd=max([max(abs(S4.alpha-S5.alpha)), ...
            max(abs(S4.Tstar-S5.Tstar)), ...
            max(abs(S4.mu-S5.mu))]);
fprintf('VIA-P5 constant-wind magnitude invariance:\n');
fprintf('   w0 = [%.1f %.1f %.1f] m/s\n',E.w0);
fprintf('   max state mismatch vs still air = %.3e\n',maxState);
fprintf('   max command mismatch            = %.3e\n',maxCmd);
assert(maxState<1e-9 && maxCmd<1e-9,'VIA-P5 failed');

%% VIA-P6 time-varying wind UUB check
S6=simulate_viA_case('varying',A.Tvarying,R,A,E);
idx=S6.t>=A.tUUB0;
linf6=max(abs(S6.e(idx,:)),[],1);
rmse6=sqrt(mean(S6.e(idx,:).^2,1));
aw=max(sqrt(sum(S6.wdot.^2,2)));

bound=[aw/A.kV, ...
       aw/(A.kGamma*R.Vmin), ...
       aw/(A.kChi*R.Vmin*cos(R.gammaMax))];

minV=min(S6.y(:,1));
maxGamma=max(abs(S6.y(:,2)));

fprintf('VIA-P6 time-varying-wind ideal-path UUB:\n');
fprintf('   max ||wdot|| = %.6f m/s^2\n',aw);
fprintf('   realized Linf [eV eGamma eChi] = [%.6f %.6e %.6e]\n',linf6);
fprintf('   conservative ISS bounds         = [%.6f %.6e %.6e]\n',bound);
fprintf('   RMSE [eV eGamma eChi]           = [%.6f %.6e %.6e]\n',rmse6);
fprintf('   min V, max |gamma|              = [%.4f m/s %.4f deg]\n', ...
        minV,rad2deg(maxGamma));
assert(all(linf6<=bound+1e-8) && minV>R.Vmin && maxGamma<R.gammaMax,'VIA-P6 failed');

%% VIA-P7 feedback-feasibility transfer under aggressive governed reference
S7=simulate_viA_case('aggressive',A.Taggressive,R,A,E);
amax=max(S7.alpha); amin=min(S7.alpha);
Tmax=max(S7.Tstar); Tmin=min(S7.Tstar);
mumax=max(abs(S7.mu));
agmin=inf;
for k=1:numel(S7.t)
    C=viA_feedback_command(S7.y(k,:).',S7.xf(k,:).',R,A);
    agmin=min(agmin,C.aGamma);
end
emax=max(abs(S7.e),[],1);
dA=max(abs(S7.alpha-S7.alphaFF));
dT=max(abs(S7.Tstar-S7.Tff));
dMu=max(abs(S7.mu-S7.muFF));

fprintf('VIA-P7 feedback-feasibility transfer (aggressive reference, still air):\n');
fprintf('   max path errors [m/s deg deg] = [%.4f %.4f %.4f]\n', ...
        emax(1),rad2deg(emax(2)),rad2deg(emax(3)));
fprintf('   alpha_c range = [%.4f %.4f] deg\n',rad2deg(amin),rad2deg(amax));
fprintf('   T* range      = [%.4f %.4f] N\n',Tmin,Tmax);
fprintf('   max |mu_c|    = %.4f deg, min a_gamma,c = %.4f m/s^2\n', ...
        rad2deg(mumax),agmin);
fprintf('   max feedback/feedforward deviations [alpha deg, T N, mu deg] = [%.4f %.4f %.4f]\n', ...
        rad2deg(dA),dT,rad2deg(dMu));

assert(amin>R.alphaLo && amax<R.alphaHi && ...
       Tmin>R.Tmin+R.DeltaT && Tmax<R.Tmax-R.DeltaT && ...
       mumax<R.muMax && agmin>R.aGammaMin,'VIA-P7 failed');

%% VIA-P8 S1 course branch crossing under ideal realization
S8=simulate_viA_case('wrap',A.Twrap,R,A,E);
maxEchi=max(abs(S8.e(:,3)));
maxJump=max(abs(diff(S8.e(:,3))));
maxRawJump=max(abs(diff(S8.cc(:,3))));
fprintf('VIA-P8 S1 course branch crossing:\n');
fprintf('   max |e_chi|          = %.4f deg\n',rad2deg(maxEchi));
fprintf('   max e_chi sample jump = %.4f deg\n',rad2deg(maxJump));
fprintf('   raw atan2 chart jump  = %.4f deg (expected near 360 at branch)\n',rad2deg(maxRawJump));
assert(maxEchi<deg2rad(1) && maxJump<deg2rad(0.2),'VIA-P8 failed');

Results=struct('RG',RG,'R',R,'A',A, ...
    'VIA1',struct('maxOrth',maxOrth,'maxForce',maxForce), ...
    'VIA2',struct('maxAlpha',maxA,'maxT',maxT,'maxMu',maxMu,'maxX',maxX,'maxN',maxN), ...
    'VIA3',struct('maxIdentity',maxId), ...
    'VIA4',S4,'VIA5',S5,'VIA6',S6,'VIA7',S7,'VIA8',S8);
end
