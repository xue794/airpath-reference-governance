function Results = test_wind_estimator_rev0()
%TEST_WIND_ESTIMATOR_REV0 Section IV wind-estimator tests WE-P1--WE-P5.
% MATLAB R2018a compatible; base MATLAB only.

P1plant=test_open_loop_rev1(); P=P1plant.Rev0.P; TR=P1plant.Rev0.trim;
E=wind_estimator_params();
fprintf('\nSection IV Wind Estimator Rev.0 validation\n');
fprintf('=========================================\n');

%% WE-P1: exact pseudo-wind reconstruction identity
maxRecon=0;
for k=1:101
    t=0.05*(k-1);
    phi=deg2rad(7*sin(0.41*t));
    theta=TR.alpha+deg2rad(3*cos(0.27*t));
    psi=deg2rad(25*sin(0.19*t));
    CIB=euler321_dcm(phi,theta,psi);
    [wI,~]=wind_estimator_test_wind(t,'varying');
    vaB=[34+2*sin(0.3*t); 1.2*sin(0.5*t); 3.0+0.6*cos(0.2*t)];
    vgB=vaB+CIB.'*wI;
    M=wind_estimator_measurement_from_truth(CIB,vgB,wI,t,E,false);
    maxRecon=max(maxRecon,norm(M.yw-wI));
end
fprintf('WE-P1 pseudo-wind reconstruction identity:\n');
fprintf('   max ||y_w-w|| = %.3e m/s\n',maxRecon);
assert(maxRecon<1e-12,'WE-P1 failed');

%% WE-P2: constant-wind convergence and analytic solution
w0=E.w0; what0=E.what0;
tg=linspace(0,E.t_constant,501); opts=odeset('RelTol',1e-11,'AbsTol',1e-12);
[tc,whatc]=ode45(@(t,x)E.Lw*(w0-x),tg,what0,opts);
errAnal=0; rel=zeros(numel(tc),1);
e0=what0-w0;
for k=1:numel(tc)
    eAnal=expm(-E.Lw*tc(k))*e0;
    eNum=whatc(k,:).'-w0;
    errAnal=max(errAnal,norm(eNum-eAnal));
    rel(k)=norm(eNum)/norm(e0);
end
idx=find(rel<=0.05,1,'first');
if isempty(idx), t95=NaN; else, t95=tc(idx); end
lmin=min(eig(E.Lw)); t95Theory=-log(0.05)/lmin;
fprintf('WE-P2 constant-wind convergence, w=[%.1f %.1f %.1f] m/s:\n',w0(1),w0(2),w0(3));
fprintf('   max analytic-solution mismatch = %.3e m/s\n',errAnal);
fprintf('   t95 numerical = %.3f s, theory = %.3f s\n',t95,t95Theory);
assert(errAnal<2e-8 && abs(t95-t95Theory)<0.03,'WE-P2 failed');

%% WE-P3: time-varying wind tracking, compare steady-state error amplitudes
ell=E.Lw(1,1); tv=linspace(0,E.t_varying,6001).';
w=zeros(numel(tv),3); for k=1:numel(tv), [wk,~]=wind_estimator_test_wind(tv(k),'varying'); w(k,:)=wk.'; end
[tt,wh]=ode45(@(t,x)E.Lw*(wind_only(t)-x),tv,E.what0,opts); %#ok<ASGLU>
e=w-wh;
fitIdx=(tv>=E.fit_window(1) & tv<=E.fit_window(2));
Omega=[0.25 0.40 0.20]; Atrue=[2.0 1.5 0.8];
Afit=zeros(3,1); Ath=zeros(3,1); ampRel=zeros(3,1);
for i=1:3
    Afit(i)=fit_sinusoid_amplitude(tv(fitIdx),e(fitIdx,i),Omega(i));
    Ath(i)=Atrue(i)*Omega(i)/sqrt(ell^2+Omega(i)^2);
    ampRel(i)=abs(Afit(i)-Ath(i))/Ath(i);
end
rmse=sqrt(mean(e(fitIdx,:).^2,1));
linf=max(abs(e(fitIdx,:)),[],1);
fprintf('WE-P3 time-varying wind tracking (steady-state fit):\n');
for i=1:3
    fprintf('   axis %d: error amplitude %.6f vs theory %.6f m/s (rel %.2f%%)\n', ...
        i,Afit(i),Ath(i),100*ampRel(i));
end
fprintf('   RMSE [N E D] = [%.4f %.4f %.4f] m/s\n',rmse(1),rmse(2),rmse(3));
fprintf('   Linf [N E D] = [%.4f %.4f %.4f] m/s\n',linf(1),linf(2),linf(3));
assert(max(ampRel)<0.025,'WE-P3 failed');

%% WE-P4: deterministic measurement-error bound
maxNw=0; maxRatio=0; maxNg=0; maxTheta=0; maxNz=0;
for k=1:1001
    t=0.03*(k-1);
    phi=deg2rad(8*sin(0.17*t)); theta=TR.alpha+deg2rad(4*cos(0.13*t)); psi=deg2rad(30*sin(0.11*t));
    CIB=euler321_dcm(phi,theta,psi); [wI,~]=wind_estimator_test_wind(t,'varying');
    vaB=[38+3*sin(0.09*t); 1.0*sin(0.21*t); 3.2+0.8*cos(0.15*t)];
    vgB=vaB+CIB.'*wI;
    M=wind_estimator_measurement_from_truth(CIB,vgB,wI,t,E,true);
    maxNw=max(maxNw,norm(M.nw)); maxRatio=max(maxRatio,norm(M.nw)/E.noise.nw_bound);
    maxNg=max(maxNg,norm(M.noise.ng)); maxTheta=max(maxTheta,norm(M.noise.dtheta)); maxNz=max(maxNz,norm(M.noise.nz));
end
fprintf('WE-P4 measurement-error robustness bound:\n');
fprintf('   max ||n_w|| = %.6f m/s\n',maxNw);
fprintf('   componentwise certified bound = %.6f m/s, max ratio = %.3f\n',E.noise.nw_bound,maxRatio);
fprintf('   generic L_T||n_z|| bound = %.6f m/s\n',E.noise.nw_bound_generic);
fprintf('   realized max ||n_g||=%.4f, theta=%.4f deg, ||n_z||=%.4f\n', ...
    maxNg,rad2deg(maxTheta),maxNz);
assert(maxNw<=E.noise.nw_bound*(1+1e-12),'WE-P4 failed');

%% WE-P5: estimator error-dynamics identity with noise and time-varying wind
maxId=0;
for k=1:501
    t=0.04*(k-1);
    phi=deg2rad(6*sin(0.23*t)); theta=TR.alpha+deg2rad(2*cos(0.29*t)); psi=deg2rad(20*sin(0.17*t));
    CIB=euler321_dcm(phi,theta,psi); [wI,wDot]=wind_estimator_test_wind(t,'varying');
    vaB=[36+2*sin(0.17*t);0.8*sin(0.31*t);3.0+0.5*cos(0.22*t)]; vgB=vaB+CIB.'*wI;
    M=wind_estimator_measurement_from_truth(CIB,vgB,wI,t,E,true);
    what=[1.5*sin(0.13*t);-0.7*cos(0.16*t);0.4*sin(0.19*t)];
    whatDot=E.Lw*(M.yw-what);
    tilde=wI-what;
    direct=wDot-whatDot;
    theory=-E.Lw*tilde+wDot-E.Lw*M.nw;
    maxId=max(maxId,norm(direct-theory));
end
fprintf('WE-P5 estimator error-dynamics identity:\n');
fprintf('   max identity residual = %.3e m/s^2\n',maxId);
assert(maxId<1e-12,'WE-P5 failed');

Results=struct('Plant',P1plant,'E',E, ...
    'P1',struct('reconstruction',maxRecon), ...
    'P2',struct('analyticMismatch',errAnal,'t95',t95,'t95Theory',t95Theory), ...
    'P3',struct('Afit',Afit,'Atheory',Ath,'ampRel',ampRel,'rmse',rmse,'linf',linf), ...
    'P4',struct('maxNw',maxNw,'bound',E.noise.nw_bound,'genericBound',E.noise.nw_bound_generic,'ratio',maxRatio), ...
    'P5',struct('identity',maxId));

    function wv=wind_only(t)
        [wv,~]=wind_estimator_test_wind(t,'varying');
    end
end
