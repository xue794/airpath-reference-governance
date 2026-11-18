function demo_section_vii()
%DEMO_SECTION_VII Diagnostic figures for Section VII.
% MATLAB R2018a compatible.

R=rg_params(); A=viA_params(); B=viB_params();
F=viC_params(); Dp=viD_params(); P7=vii_params(); Ew=wind_estimator_params();

S=simulate_section_vii_case('aggressive',P7.Taggressive,R,A,F,Dp,P7,Ew);

figure;
subplot(3,1,1);
plot(S.t,rad2deg(S.deltaCmd(:,1)),'k--','LineWidth',1.0); hold on;
plot(S.t,rad2deg(S.deltaAct(:,1)),'LineWidth',1.2); grid on;
ylabel('delta_a [deg]'); legend('cmd','act');
subplot(3,1,2);
plot(S.t,rad2deg(S.deltaCmd(:,2)),'k--','LineWidth',1.0); hold on;
plot(S.t,rad2deg(S.deltaAct(:,2)),'LineWidth',1.2); grid on;
ylabel('delta_e [deg]'); legend('cmd','act');
subplot(3,1,3);
plot(S.t,rad2deg(S.deltaCmd(:,3)),'k--','LineWidth',1.0); hold on;
plot(S.t,rad2deg(S.deltaAct(:,3)),'LineWidth',1.2); grid on;
ylabel('delta_r [deg]'); xlabel('t [s]'); legend('cmd','act');
title('Physical control-surface interface');

figure;
plot(S.t,sqrt(sum(S.dAlloc.^2,2)),'LineWidth',1.1); hold on;
plot(S.t,sqrt(sum(S.dServo.^2,2)),'LineWidth',1.1);
plot(S.t,sqrt(sum(S.dTau.^2,2)),'LineWidth',1.2);
plot(S.t,P7.dTauBudget*ones(size(S.t)),'k--');
grid on; xlabel('t [s]'); ylabel('moment residual [N m]');
legend('|d_{alloc}|','|d_{servo}|','|d_tau|','VI-D steady residual budget');
title('Allocation and actuator moment residuals');

figure;
plot(S.t,rad2deg(S.thetaF),'LineWidth',1.0); hold on;
plot(S.t,rad2deg(S.thetaR),'LineWidth',1.0);
plot(S.t,rad2deg(S.thetaF+S.thetaR),'LineWidth',1.2);
plot(S.t,rad2deg(B.thetaTube)*ones(size(S.t)),'k--');
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('theta_f','theta_R','total','certified tube');
title('Attitude certificate after physical surface servos');

figure;
subplot(2,1,1);
plot(S.t,S.Tcmd,'k--','LineWidth',1.0); hold on;
plot(S.t,S.Tact,'LineWidth',1.2); grid on;
ylabel('T [N]'); legend('T_cmd','T_act');
subplot(2,1,2);
plot(S.t,sqrt(sum(S.fDelta.^2,2)),'LineWidth',1.2); grid on;
ylabel('|d_F,m,surf| [N]'); xlabel('t [s]');
title('Propulsor interface and omitted surface-force residual');

figure;
subplot(2,1,1);
plot(S.t,rad2deg(S.alpha),'LineWidth',1.2); grid on; ylabel('alpha [deg]');
subplot(2,1,2);
plot(S.t,rad2deg(S.beta),'LineWidth',1.2); grid on;
ylabel('beta [deg]'); xlabel('t [s]');
title('Airdata under dynamically allocated surfaces');
end
