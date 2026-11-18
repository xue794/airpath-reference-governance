function Results = test_reference_governor_rev0()
%TEST_REFERENCE_GOVERNOR_REV0 Section V module tests RG-P1--RG-P8.

SMC = test_smc_guidance_rev0();
R=rg_params();

fprintf('\nSection V Certified Constrained Air-Path Reference Governor Rev.0\n');
fprintf('================================================================\n');

%% RG-P1 raw-command chart and S1 inverse consistency
maxInv=0; maxAngle=0;
for k=1:801
    t=0.071*(k-1);
    [vc,c0]=rg_raw_test_command(t,'wrap',R);
    [cc,info]=rg_raw_command_chart(vc,R);
    vr=rg_airpath_to_velocity(cc);
    maxInv=max(maxInv,norm(vr-vc));
    maxAngle=max(maxAngle,abs(rg_angle_error(cc(3),c0(3))));
    if ~info.regular, error('RG-P1 chart regularity lost'); end
end
fprintf('RG-P1 raw-command chart / S1 consistency:\n');
fprintf('   max velocity reconstruction error = %.3e m/s\n',maxInv);
fprintf('   max wrapped course error          = %.3e rad\n',maxAngle);
assert(maxInv<1e-12 && maxAngle<1e-12,'RG-P1 failed');

%% RG-P2 monotonicity scan and backup feedforward root
minNa=inf;
for V=linspace(R.Vmin,R.Vmax,17)
    for gam=linspace(-R.gammaMax,R.gammaMax,9)
        for nuV=linspace(-R.nuMax(1),R.nuMax(1),5)
            X=R.P.m*(nuV+R.g*sin(gam));
            for a=linspace(R.alphaLo,R.alphaHi,37)
                M=rg_N_map(a,V,X,R);
                minNa=min(minNa,M.alpha);
            end
        end
    end
end
xb=[R.Vb;R.gammaB;0;0;0;0];
Fb=rg_feedforward(xb,R);
fprintf('RG-P2 combined-force inversion:\n');
fprintf('   grid min dN/dalpha = %.3f N/rad\n',minNa);
fprintf('   backup alpha_ff    = %.6f deg\n',rad2deg(Fb.alpha));
fprintf('   backup T_ff        = %.6f N\n',Fb.T);
fprintf('   root residual      = %.3e N\n',abs(Fb.rootInfo.residual));
assert(minNa>0 && Fb.rootInfo.bracket && abs(Fb.rootInfo.residual)<1e-10,'RG-P2 failed');

%% RG-P3 IFT feedforward-rate and affine-constraint derivative audit
rng(3,'twister');
maxAlphaRateErr=0; maxTRateErr=0; maxAffine=0;
for k=1:200
    dx=[0.30*(2*rand-1);
        deg2rad(0.5)*(2*rand-1);
        0;
        0.15*(2*rand-1);
        deg2rad(0.5)*(2*rand-1);
        deg2rad(1.0)*(2*rand-1)];
    x=xb+dx;
    u=[0.25*(2*rand-1);deg2rad(1.5)*(2*rand-1);deg2rad(2.0)*(2*rand-1)];
    F0=rg_feedforward(x,R);
    D0=rg_feedforward_rate(x,u,R);
    dt=1e-6;
    xdot=[x(4:6);u];
    F1=rg_feedforward(x+dt*xdot,R);
    adFD=rg_angle_error(F1.alpha,F0.alpha)/dt;
    TdFD=(F1.T-F0.T)/dt;
    maxAlphaRateErr=max(maxAlphaRateErr,abs(adFD-D0.alphaDot));
    maxTRateErr=max(maxTRateErr,abs(TdFD-D0.Tdot));

    C=rg_constraint_affine(x,R);
    phi1=rg_constraint_direct(x,u,R);
    phi2=C.c+C.G*u;
    maxAffine=max(maxAffine,max(abs(phi1-phi2)));
end
fprintf('RG-P3 feedforward derivative / affine constraint audit:\n');
fprintf('   max alpha_dot finite-difference error = %.3e rad/s\n',maxAlphaRateErr);
fprintf('   max T_dot finite-difference error     = %.3e N/s\n',maxTRateErr);
fprintf('   max affine-evaluation residual        = %.3e\n',maxAffine);
assert(maxAlphaRateErr<2e-5 && maxTRateErr<5e-4 && maxAffine<1e-12,'RG-P3 failed');

%% RG-P4 backup trim, Lyapunov equation, and dense local K_B audit
lyapRes=norm(R.AB'*R.PB+R.PB*R.AB+R.QB,'fro');
Cb=rg_constraint_affine(xb,R);
sigb=Cb.c+Cb.G*rg_backup_control(xb,R);
Audit=rg_backup_sample_audit(R,6000);
fprintf('RG-P4 strict backup set audit:\n');
fprintf('   Lyapunov-equation residual = %.3e\n',lyapRes);
fprintf('   trim min backup slack      = %.6f\n',min(sigb));
fprintf('   sampled K_B min slack      = %.6f (%s)\n',Audit.minSigma,Audit.worstName);
fprintf('   sampled min h, psi1        = [%.6f %.6f]\n',Audit.minH,Audit.minPsi);
fprintf('   sampled min N-margin       = %.6f N\n',Audit.minRootMargin);
assert(lyapRes<1e-10 && min(sigb)>0 && Audit.minSigma>0 && ...
       Audit.minH>=0 && Audit.minPsi>=0 && Audit.minRootMargin>0,'RG-P4 failed');

%% RG-P5 pointwise scalar-governor feasibility
rng(5,'twister');
[Vd,Dd]=eig(R.PB); PinvHalf=Vd*diag(1./sqrt(diag(Dd)))*Vd';
minPhi=inf; minLam=1; maxLam=0;
for k=1:2000
    d=randn(5,1); d=d/norm(d);
    xi=PinvHalf*(sqrt(0.95*R.rhoB)*rand()^(1/5)*d);
    x=[R.Vb+xi(1);R.gammaB+xi(2);0;xi(3);xi(4);xi(5)];
    u0=rg_backup_control(x,R)+[4*(2*rand-1);deg2rad(30)*(2*rand-1);deg2rad(40)*(2*rand-1)];
    O=rg_scalar_governor(x,u0,R);
    minPhi=min(minPhi,min(O.phi));
    minLam=min(minLam,O.lambda); maxLam=max(maxLam,O.lambda);
end
fprintf('RG-P5 scalar-governor pointwise feasibility:\n');
fprintf('   minimum realized phi_j = %.3e\n',minPhi);
fprintf('   lambda range           = [%.6f %.6f]\n',minLam,maxLam);
assert(minPhi>-1e-10 && minLam>0 && maxLam<=1+1e-12,'RG-P5 failed');

%% RG-P6 exact pass-through when u0 is feasible
S6=simulate_reference_governor('small',R.Tsmall,R);
minLam6=min(S6.lambda); maxDiff=0; minPhi6=inf;
for k=1:numel(S6.t)
    x=S6.x(k,:).';
    [vc,~]=rg_raw_test_command(S6.t(k),'small',R);
    [cc,~]=rg_raw_command_chart(vc,R);
    u0=rg_nominal_control(x,cc,R);
    O=rg_scalar_governor(x,u0,R);
    maxDiff=max(maxDiff,norm(O.u-u0));
    minPhi6=min(minPhi6,min(O.phi));
end
fprintf('RG-P6 feasible-command pass-through:\n');
fprintf('   min lambda        = %.12f\n',minLam6);
fprintf('   max ||u*-u0||     = %.3e\n',maxDiff);
fprintf('   min realized phi  = %.6f\n',minPhi6);
assert(abs(minLam6-1)<1e-12 && maxDiff<1e-10 && minPhi6>0,'RG-P6 failed');

%% RG-P7 aggressive command: recursive feasibility and active governor
S7=simulate_reference_governor('aggressive',R.Taggressive,R);
minPhi7=min(S7.minphi);
minH7=min(S7.hB);
minLam7=min(S7.lambda); maxLam7=max(S7.lambda);
maxJump=max(abs(diff(S7.lambda)));
amin7=min(S7.alpha); amax7=max(S7.alpha);
Tmin7=min(S7.Tff); Tmax7=max(S7.Tff);
mumax7=max(abs(S7.mu));
fprintf('RG-P7 aggressive raw-command simulation:\n');
fprintf('   lambda range         = [%.6f %.6f], max sample jump %.4f\n',minLam7,maxLam7,maxJump);
fprintf('   min phi_j, min h_B   = [%.3e %.3e]\n',minPhi7,minH7);
fprintf('   alpha_ff range       = [%.4f %.4f] deg\n',rad2deg(amin7),rad2deg(amax7));
fprintf('   T_ff range           = [%.4f %.4f] N\n',Tmin7,Tmax7);
fprintf('   max |mu_ff|          = %.4f deg\n',rad2deg(mumax7));
assert(minPhi7>-1e-8 && minH7>-1e-8 && minLam7>0 && minLam7<0.95 && ...
       amax7<R.alphaHi && amin7>R.alphaLo && ...
       Tmin7>R.Tmin+R.DeltaT && Tmax7<R.Tmax-R.DeltaT && mumax7<R.muMax, ...
       'RG-P7 failed');

%% RG-P8 course branch crossing
S8=simulate_reference_governor('wrap',R.Twrap,R);
maxEtaJump=max(abs(diff(S8.etaChi)));
maxEta=max(abs(S8.etaChi));
minPhi8=min(S8.minphi);
fprintf('RG-P8 S1 course branch crossing:\n');
fprintf('   max wrapped course error = %.4f deg\n',rad2deg(maxEta));
fprintf('   max eta_chi sample jump  = %.4f deg\n',rad2deg(maxEtaJump));
fprintf('   min realized phi_j       = %.3e\n',minPhi8);
assert(maxEta<pi-R.course_chart_margin && maxEtaJump<deg2rad(2) && minPhi8>-1e-8,'RG-P8 failed');

Results=struct('SMC',SMC,'R',R,'RG1',struct('maxInv',maxInv,'maxAngle',maxAngle), ...
 'RG2',struct('minNalpha',minNa,'backupFeedforward',Fb), ...
 'RG3',struct('alphaRateErr',maxAlphaRateErr,'TRateErr',maxTRateErr,'affineErr',maxAffine), ...
 'RG4',struct('lyapResidual',lyapRes,'trimSlack',min(sigb),'audit',Audit), ...
 'RG5',struct('minPhi',minPhi,'minLambda',minLam,'maxLambda',maxLam), ...
 'RG6',struct('sim',S6,'minLambda',minLam6,'maxPassThroughError',maxDiff), ...
 'RG7',S7,'RG8',S8);
end
