function demo_section_viC()
%DEMO_SECTION_VIC Diagnostic figures for Section VI-C.
% MATLAB R2018a compatible; no tiledlayout/yline required.

R=rg_params(); A=viA_params(); F=viC_params(); Ew=wind_estimator_params();

S=simulate_viC_case('aggressive',A.Taggressive,R,A,F,Ew);

figure;
subplot(2,1,1);
plot(S.t,rad2deg(S.theta),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(F.thetaFMax)*ones(size(S.t)),'k--');
grid on; ylabel('theta_f [deg]');
legend('theta_f','certified 0.30-deg tube');
title('SO(3) command-filter attitude error');
subplot(2,1,2);
plot(S.t,S.Psi,'LineWidth',1.2); hold on;
plot(S.t,F.PsiFMax*ones(size(S.t)),'k--');
grid on; ylabel('Psi_f'); xlabel('t [s]');
legend('Psi_f','tube-equivalent bound');

figure;
plot(S.t,sqrt(sum(S.omegaC.^2,2)),'k--','LineWidth',1.0); hold on;
plot(S.t,sqrt(sum(S.omegaR.^2,2)),'LineWidth',1.2);
plot(S.t,sqrt(sum(S.xi.^2,2)),'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('angular rate [rad/s]');
legend('|omega_c|','|omega_r|','|xi_f|');
title('Command-filter angular-rate variables');

figure;
subplot(3,1,1);
plot(S.t,S.e(:,1),'LineWidth',1.2); grid on; ylabel('e_{Rf,1}');
subplot(3,1,2);
plot(S.t,S.e(:,2),'LineWidth',1.2); grid on; ylabel('e_{Rf,2}');
subplot(3,1,3);
plot(S.t,S.e(:,3),'LineWidth',1.2); grid on; ylabel('e_{Rf,3}'); xlabel('t [s]');

% Persistent rotating-command demonstration.
tg=linspace(0,F.spinT,4001);
x0=[1;0;0;0;0;0;0];
opts=odeset('RelTol',F.RelTol,'AbsTol',F.AbsTol,'MaxStep',0.01);
[t,x]=ode45(@(tt,xx)viC_spin_rhs(tt,xx,F),tg,x0,opts);
theta=zeros(numel(t),1);
for k=1:numel(t)
    theta(k)=dcm_distance(viB_Rz(F.spinRate*t(k)),quat_to_dcm(x(k,1:4).'));
end
pred=asin(F.KomegaF(3,3)*F.spinRate/F.kRf);

figure;
plot(t,rad2deg(theta),'LineWidth',1.2); hold on;
plot(t,rad2deg(pred)*ones(size(t)),'k--');
grid on; xlabel('t [s]'); ylabel('theta_f [deg]');
legend('moving-command filter error','exact steady lag');
title('Persistent moving reference produces bounded nonzero filter error');

% S1 branch-crossing demonstration.
W=simulate_viC_case('wrap',A.Twrap,R,A,F,Ew);
figure;
subplot(2,1,1);
plot(W.t,rad2deg(W.Up.cc(:,3)),'k--','LineWidth',1.0); hold on;
plot(W.t,rad2deg(W.Up.y(:,3)),'LineWidth',1.2); grid on;
ylabel('course [deg]');
legend('raw atan2 chi_c','continuous air-path chi');
title('S1 branch crossing');
subplot(2,1,2);
plot(W.t,rad2deg(W.theta),'LineWidth',1.2); hold on;
plot(W.t,rad2deg(F.thetaFMax)*ones(size(W.t)),'k--');
grid on; ylabel('theta_f [deg]'); xlabel('t [s]');
legend('theta_f','certified tube');
end
