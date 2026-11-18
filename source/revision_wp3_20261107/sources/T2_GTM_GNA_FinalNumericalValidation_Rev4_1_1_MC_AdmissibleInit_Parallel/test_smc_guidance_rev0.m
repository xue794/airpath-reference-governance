function Results = test_smc_guidance_rev0()
%TEST_SMC_GUIDANCE_REV0 Section IV guidance module tests SMC-P1--SMC-P7.
% MATLAB R2018a compatible; base MATLAB only.

WE = test_wind_estimator_rev0();
E = WE.E;
G = smc_guidance_params();

fprintf('\nSection IV Wind-Compensated Integral-Error SMC Guidance Rev.0\n');
fprintf('============================================================\n');

%% SMC-P1: exact sliding-dynamics identity
maxId=0;
for k=1:401
    t=0.037*(k-1);
    [rd,~]=smc_guidance_reference(t);
    e=[5*sin(0.13*t);-4*cos(0.17*t);3*sin(0.19*t+0.2)];
    ze=[1.1*cos(0.11*t);-0.7*sin(0.14*t);0.4*cos(0.09*t)];
    what=[1.5*sin(0.10*t);-0.8*cos(0.12*t);0.3*sin(0.08*t)];
    z=[rd+e;ze;what];
    O=smc_guidance_command(t,z(1:3),ze,what,G);
    S=smc_guidance_test_signals(t,'bound',E);
    zDot=smc_guidance_augmented_rhs(t,z,G,E,'bound');
    erDot=zDot(1:3)-O.rdDot;
    zeDot=zDot(4:6);
    sDotDirect=erDot+G.Lambda*zeDot;
    sDotTheory=-G.Ks*O.s-G.Kc.*smc_sat(O.s./G.eps)+(S.w-what)+S.ev;
    maxId=max(maxId,norm(sDotDirect-sDotTheory));
end
fprintf('SMC-P1 exact sliding-dynamics identity:\n');
fprintf('   max identity residual = %.3e m/s\n',maxId);
assert(maxId<1e-12,'SMC-P1 failed');

%% SMC-P2: still-air exact tracking
R2=simulate_smc_guidance_case('still',G.t_still,G,E);
idx=R2.t>=50;
lateE=max(abs(R2.er(idx,:)),[],1);
lateS=max(abs(R2.s(idx,:)),[],1);
fprintf('SMC-P2 still-air exact tracking:\n');
fprintf('   late Linf |e_r| [N E D] = [%.3e %.3e %.3e] m\n',lateE);
fprintf('   late Linf |s|   [N E D] = [%.3e %.3e %.3e] m\n',lateS);
assert(max(lateE)<1e-6 && max(lateS)<1e-8,'SMC-P2 failed');

%% SMC-P3: constant-wind compensation
R3=simulate_smc_guidance_case('constant',G.t_constant,G,E);
idx=R3.t>=50;
lateE3=max(abs(R3.er(idx,:)),[],1);
lateS3=max(abs(R3.s(idx,:)),[],1);
lateW3=max(abs(R3.tildew(idx,:)),[],1);
fprintf('SMC-P3 constant-wind compensation:\n');
fprintf('   late Linf |e_r|    = [%.3e %.3e %.3e] m\n',lateE3);
fprintf('   late Linf |s|      = [%.3e %.3e %.3e] m\n',lateS3);
fprintf('   late Linf |w-what| = [%.3e %.3e %.3e] m/s\n',lateW3);
assert(max(lateE3)<1e-6 && max(lateS3)<1e-8 && max(lateW3)<1e-8,'SMC-P3 failed');

%% SMC-P4: time-varying wind UUB / exact sinusoidal amplitude check
R4=simulate_smc_guidance_case('varying',G.t_varying,G,E);
idx=R4.t>=G.fit_window(1) & R4.t<=G.fit_window(2);
Omega=[0.25 0.40 0.20]; Aw=[2.0 1.5 0.8];
ell=E.Lw(1,1);
kEff=diag(G.Ks)+G.Kc./G.eps;
Afit=zeros(3,1); Ath=zeros(3,1); rel=zeros(3,1);
for i=1:3
    Afit(i)=fit_sinusoid_amplitude(R4.t(idx),R4.s(idx,i),Omega(i));
    Ae=Aw(i)*Omega(i)/sqrt(ell^2+Omega(i)^2);
    Ath(i)=Ae/sqrt(kEff(i)^2+Omega(i)^2);
    rel(i)=abs(Afit(i)-Ath(i))/Ath(i);
end
rmseE=sqrt(mean(R4.er(idx,:).^2,1));
linfS=max(abs(R4.s(idx,:)),[],1);
fprintf('SMC-P4 time-varying-wind UUB (steady-state fit):\n');
for i=1:3
    fprintf('   axis %d: |s| amplitude %.6f vs theory %.6f m (rel %.2f%%)\n', ...
        i,Afit(i),Ath(i),100*rel(i));
end
fprintf('   RMSE e_r [N E D] = [%.4f %.4f %.4f] m\n',rmseE);
fprintf('   Linf s   [N E D] = [%.4f %.4f %.4f] m\n',linfS);
assert(max(rel)<0.025 && all(linfS<G.eps.'),'SMC-P4 failed');

%% SMC-P5: governor-mismatch cancellation
R50=simulate_smc_guidance_case('varying',G.t_mismatch,G,E);
R5g=simulate_smc_guidance_case('mismatch',G.t_mismatch,G,E);
sMismatch=max(sqrt(sum((R50.s-R5g.s).^2,2)));
eDifference=max(sqrt(sum((R50.er-R5g.er).^2,2)));
fprintf('SMC-P5 governor-mismatch cancellation:\n');
fprintf('   max ||s(e_g)-s(0)|| = %.3e m\n',sMismatch);
fprintf('   max ||e_r(e_g)-e_r(0)|| = %.6f m (expected nonzero)\n',eDifference);
assert(sMismatch<1e-4 && eDifference>0.1,'SMC-P5 failed');

%% SMC-P6: componentwise position ultimate bound

% MATLAB R2018a note:
% In SMC-P6, rs is a 3x1 column vector. Compare it with G.eps (also 3x1).
% Using G.eps.' would invoke implicit expansion to a 3x3 logical matrix,
% so all(...) would return a 1x3 vector and the && operator would fail.

R6=simulate_smc_guidance_case('bound',G.t_bound,G,E);
idx=R6.t>=G.bound_window(1) & R6.t<=G.bound_window(2);
rs=max(abs(R6.s(idx,:)),[],1).';
rg=max(abs(R6.eg(idx,:)),[],1).';
re=max(abs(R6.er(idx,:)),[],1).';
bound=2*rs+rg./diag(G.Lambda);
dbar=max(abs(R6.tildew(idx,:)+R6.ev(idx,:)),[],1).';
refinedS=dbar./kEff;
ratio=re./bound;
fprintf('SMC-P6 boundary-layer position bound:\n');
fprintf('   realized Linf e_r = [%.4f %.4f %.4f] m\n',re);
fprintf('   theorem bound     = [%.4f %.4f %.4f] m\n',bound);
fprintf('   realized/bound    = [%.3f %.3f %.3f]\n',ratio);
fprintf('   Linf s            = [%.4f %.4f %.4f] m\n',rs);
fprintf('   dbar/k_eff        = [%.4f %.4f %.4f] m\n',refinedS);
assert(all(re<=bound+1e-8) && all(rs<G.eps) && all(rs<=1.05*refinedS),'SMC-P6 failed');

%% SMC-P7: estimator-measurement-error robustness
R7=simulate_smc_guidance_case('noisy',G.t_bound,G,E);
idx=R7.t>=G.bound_window(1) & R7.t<=G.bound_window(2);
rmse7=sqrt(mean(R7.er(idx,:).^2,1));
linf7=max(abs(R7.er(idx,:)),[],1);
linfTw7=max(abs(R7.tildew(idx,:)),[],1);
linfS7=max(abs(R7.s(idx,:)),[],1);
nwNorm=sqrt(sum(R7.nw(idx,:).^2,2));
maxNw=max(nwNorm);
fprintf('SMC-P7 bounded pseudo-wind measurement-error robustness:\n');
fprintf('   max ||n_w|| = %.4f m/s\n',maxNw);
fprintf('   RMSE e_r [N E D] = [%.4f %.4f %.4f] m\n',rmse7);
fprintf('   Linf e_r [N E D] = [%.4f %.4f %.4f] m\n',linf7);
fprintf('   Linf tilde-w      = [%.4f %.4f %.4f] m/s\n',linfTw7);
fprintf('   Linf s            = [%.4f %.4f %.4f] m\n',linfS7);
assert(maxNw<E.noise.nw_bound && all(linfS7<G.eps.') && max(linf7)<2.0,'SMC-P7 failed');

Results=struct('WindEstimator',WE,'G',G, ...
    'P1',struct('identity',maxId), ...
    'P2',struct('lateE',lateE,'lateS',lateS), ...
    'P3',struct('lateE',lateE3,'lateS',lateS3,'lateW',lateW3), ...
    'P4',struct('Afit',Afit,'Atheory',Ath,'rel',rel,'rmseE',rmseE,'linfS',linfS), ...
    'P5',struct('sMismatch',sMismatch,'eDifference',eDifference), ...
    'P6',struct('rs',rs,'rg',rg,'realizedE',re,'bound',bound,'ratio',ratio,'refinedS',refinedS), ...
    'P7',struct('maxNw',maxNw,'rmseE',rmse7,'linfE',linf7,'linfTw',linfTw7,'linfS',linfS7));
end
