function demo_section_viD()
%DEMO_SECTION_VID Diagnostic figures for Section VI-D.
% MATLAB R2018a compatible.

R=rg_params(); A=viA_params(); B=viB_params();
F=viC_params(); Dp=viD_params(); Ew=wind_estimator_params();

S=simulate_viD_case('aggressive_disturbed',Dp.Taggressive,R,A,F,Dp,Ew);

figure;
plot(S.t,rad2deg(S.thetaF),'LineWidth',1.1); hold on;
plot(S.t,rad2deg(S.thetaR),'LineWidth',1.1);
plot(S.t,rad2deg(S.thetaF+S.thetaR),'LineWidth',1.2);
plot(S.t,rad2deg(B.thetaTube)*ones(size(S.t)),'k--');
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('theta_f','theta_R','theta_f+theta_R','certified 0.70-deg tube');
title('Dynamic realization of the complete VI-B attitude tube');

figure;
subplot(2,1,1);
plot(S.t,S.eR,'LineWidth',1.1); grid on;
ylabel('e_R'); legend('1','2','3');
title('Rotational attitude error');
subplot(2,1,2);
plot(S.t,S.eOmega,'LineWidth',1.1); grid on;
ylabel('e_omega [rad/s]'); xlabel('t [s]'); legend('1','2','3');

figure;
subplot(2,1,1);
plot(S.t,S.tauC,'LineWidth',1.1); grid on;
ylabel('tau_c [N m]'); legend('roll','pitch','yaw');
title('Ideal PHS virtual torque');
subplot(2,1,2);
plot(S.t,S.dTau,'LineWidth',1.1); grid on;
ylabel('d_tau [N m]'); xlabel('t [s]'); legend('roll','pitch','yaw');

figure;
subplot(3,1,1);
plot(S.t,rad2deg(S.alpha),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(R.alphaPhysicalMin)*ones(size(S.t)),'k--');
plot(S.t,rad2deg(R.alphaPhysicalMax)*ones(size(S.t)),'k--');
grid on; ylabel('alpha [deg]');
subplot(3,1,2);
plot(S.t,rad2deg(S.beta),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(B.betaBar)*ones(size(S.t)),'k--');
plot(S.t,-rad2deg(B.betaBar)*ones(size(S.t)),'k--');
grid on; ylabel('beta [deg]');
subplot(3,1,3);
plot(S.t,S.Tphys,'LineWidth',1.2); hold on;
plot(S.t,(R.Tmin+R.DeltaT)*ones(size(S.t)),'k--');
plot(S.t,(R.Tmax-R.DeltaT)*ones(size(S.t)),'k--');
grid on; ylabel('T [N]'); xlabel('t [s]');
title('Physical geometry under dynamic Q_f and Q_R');

W=simulate_viD_case('wrap',Dp.Twrap,R,A,F,Dp,Ew);
figure;
subplot(2,1,1);
plot(W.t,rad2deg(W.Up.cc(:,3)),'k--','LineWidth',1.0); hold on;
plot(W.t,rad2deg(W.Up.y(:,3)),'LineWidth',1.2);
grid on; ylabel('course [deg]');
legend('raw atan2 chi_c','continuous air-path chi');
title('S1 branch crossing through VI-D');
subplot(2,1,2);
plot(W.t,rad2deg(W.thetaR),'LineWidth',1.2); hold on;
plot(W.t,rad2deg(Dp.thetaRMax)*ones(size(W.t)),'k--');
grid on; ylabel('theta_R [deg]'); xlabel('t [s]');
legend('theta_R','certified 0.40-deg tube');
end
