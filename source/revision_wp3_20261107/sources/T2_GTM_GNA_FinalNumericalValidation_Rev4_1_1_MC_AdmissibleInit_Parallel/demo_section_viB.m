function demo_section_viB()
%DEMO_SECTION_VIB Section VI-B diagnostic plots; MATLAB R2018a compatible.

R=rg_params(); A=viA_params(); B=viB_params(); E=wind_estimator_params();
S=simulate_viA_case('aggressive',A.Taggressive,R,A,E);

n=numel(S.t);
alpha=zeros(n,1); beta=zeros(n,1); Tphys=zeros(n,1);
thetaF=zeros(n,1); thetaR=zeros(n,1); kappa=zeros(n,1);

for k=1:n
    y=S.y(k,:).';
    Syn=viB_synthetic_relative_rotations(S.t(k),B);
    G=viB_geometry(y(2),y(3),S.mu(k),S.alpha(k),Syn.Qf,Syn.QR);
    alpha(k)=G.alpha; beta(k)=G.beta; kappa(k)=G.kappa;
    Tphys(k)=S.Fparallel(k)/G.kappa;
    thetaF(k)=G.thetaF; thetaR(k)=G.thetaR;
end

figure;
subplot(2,1,1);
plot(S.t,rad2deg(thetaF),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(thetaR),'LineWidth',1.2);
plot(S.t,rad2deg(thetaF+thetaR),'k--','LineWidth',1.0);
plot(S.t,rad2deg(B.thetaTube)*ones(size(S.t)),'k:');
grid on; ylabel('angle [deg]');
legend('theta_f','theta_R','theta_f+theta_R','certified tube');
title('Synthetic relative-attitude errors');
subplot(2,1,2);
plot(S.t,kappa,'LineWidth',1.2); hold on;
plot(S.t,B.kappaMin*ones(size(S.t)),'k--');
grid on; ylabel('kappa_T'); xlabel('t [s]');
legend('kappa_T','certified lower bound');

figure;
subplot(2,1,1);
plot(S.t,rad2deg(S.alpha),'k--','LineWidth',1.0); hold on;
plot(S.t,rad2deg(alpha),'LineWidth',1.2);
plot(S.t,rad2deg(R.alphaPhysicalMax)*ones(size(S.t)),'k:');
plot(S.t,rad2deg(R.alphaPhysicalMin)*ones(size(S.t)),'k:');
grid on; ylabel('alpha [deg]');
legend('alpha_c','actual alpha','physical bounds');
subplot(2,1,2);
plot(S.t,rad2deg(beta),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(B.betaBar)*ones(size(S.t)),'k--');
plot(S.t,-rad2deg(B.betaBar)*ones(size(S.t)),'k--');
grid on; ylabel('beta [deg]'); xlabel('t [s]');
legend('beta','+beta bound','-beta bound');

figure;
plot(S.t,S.Tstar,'k--','LineWidth',1.0); hold on;
plot(S.t,Tphys,'LineWidth',1.2);
plot(S.t,(R.Tmin+R.DeltaT)*ones(size(S.t)),'k:');
plot(S.t,(R.Tmax-R.DeltaT)*ones(size(S.t)),'k:');
grid on; xlabel('t [s]'); ylabel('thrust [N]');
legend('ideal T* (alpha=alpha_c,beta=0)','geometric T*','reference interior bounds');
title('Body-axis thrust geometry transfer');

% S1 branch representation test.
W=simulate_viA_case('wrap',A.Twrap,R,A,E);
step=zeros(numel(W.t)-1,1);
Cprev=[];
for k=1:numel(W.t)
    y=W.y(k,:).';
    Cc=viB_command_attitude(y(2),y(3),W.mu(k),W.alpha(k));
    if ~isempty(Cprev)
        step(k-1)=dcm_distance(Cprev,Cc);
    end
    Cprev=Cc;
end
figure;
subplot(2,1,1);
plot(W.t,rad2deg(W.cc(:,3)),'k--','LineWidth',1.0); hold on;
plot(W.t,rad2deg(W.y(:,3)),'LineWidth',1.2); grid on;
ylabel('course [deg]'); legend('raw chi_c from atan2','actual chi');
title('S1 branch crossing');
subplot(2,1,2);
plot(W.t(2:end),rad2deg(step),'LineWidth',1.2); grid on;
ylabel('DCM step [deg]'); xlabel('t [s]');
title('Command-attitude geodesic step remains continuous');
end
