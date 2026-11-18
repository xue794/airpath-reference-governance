function demo_complete_6dof_rev0(Results)
%DEMO_COMPLETE_6DOF_REV0 Diagnostic plots for the fully integrated loop.
% Pass Results from test_complete_6dof_rev0 to avoid rerunning long cases.

if nargin<1 || isempty(Results)
    C=cl6_params();
    V=simulate_complete_6dof_case('varying',C);
    N=simulate_complete_6dof_case('noisy',C);
else
    C=Results.C;
    V=Results.varying;
    N=Results.noisy;
end

figure;
subplot(3,1,1);
plot(V.t,V.Y.er(:,1),'LineWidth',1.1); grid on; ylabel('e_N [m]');
subplot(3,1,2);
plot(V.t,V.Y.er(:,2),'LineWidth',1.1); grid on; ylabel('e_E [m]');
subplot(3,1,3);
plot(V.t,V.Y.er(:,3),'LineWidth',1.1); grid on; ylabel('e_D [m]'); xlabel('t [s]');
title('Complete 6DOF ground-position tracking errors');

figure;
plot(V.t,rad2deg(V.Y.thetaF),'LineWidth',1.1); hold on;
plot(V.t,rad2deg(V.Y.thetaR),'LineWidth',1.1);
plot(V.t,rad2deg(V.Y.thetaF+V.Y.thetaR),'LineWidth',1.2);
plot(V.t,rad2deg(C.B.thetaTube)*ones(size(V.t)),'k--');
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('theta_f','theta_R','total','0.70-deg tube');
title('Complete dynamic attitude certificate');

figure;
subplot(3,1,1);
plot(V.t,V.Y.V,'LineWidth',1.1); hold on;
plot(V.t,C.R.Vmin*ones(size(V.t)),'k--');
plot(V.t,C.R.Vmax*ones(size(V.t)),'k--');
grid on; ylabel('V [m/s]');
subplot(3,1,2);
plot(V.t,rad2deg(V.Y.alpha),'LineWidth',1.1); grid on; ylabel('alpha [deg]');
subplot(3,1,3);
plot(V.t,rad2deg(V.Y.beta),'LineWidth',1.1); grid on;
ylabel('beta [deg]'); xlabel('t [s]');
title('Actual air-relative states from the nonlinear 6DOF plant');

figure;
subplot(3,1,1);
plot(V.t,rad2deg(V.Y.deltaCmd(:,1)),'k--'); hold on;
plot(V.t,rad2deg(V.Y.deltaAct(:,1)),'LineWidth',1.1); grid on;
ylabel('delta_a [deg]'); legend('cmd','act');
subplot(3,1,2);
plot(V.t,rad2deg(V.Y.deltaCmd(:,2)),'k--'); hold on;
plot(V.t,rad2deg(V.Y.deltaAct(:,2)),'LineWidth',1.1); grid on;
ylabel('delta_e [deg]'); legend('cmd','act');
subplot(3,1,3);
plot(V.t,rad2deg(V.Y.deltaCmd(:,3)),'k--'); hold on;
plot(V.t,rad2deg(V.Y.deltaAct(:,3)),'LineWidth',1.1); grid on;
ylabel('delta_r [deg]'); xlabel('t [s]'); legend('cmd','act');
title('Physical surface commands and servo states');

figure;
plot(V.t,sqrt(sum(V.Y.dTauPhysical.^2,2)),'LineWidth',1.1); hold on;
plot(V.t,C.P7.dTauBudget*ones(size(V.t)),'k--');
grid on; xlabel('t [s]'); ylabel('|d_tau,physical| [N m]');
legend('physical rotational residual','old VI-D scalar budget');
title('Full-system rotational residual: actuator state is retained explicitly');

figure;
subplot(2,1,1);
plot(V.t,V.Y.Tcmd,'k--'); hold on;
plot(V.t,V.Y.Tact,'LineWidth',1.1); grid on;
ylabel('T [N]'); legend('cmd','act');
subplot(2,1,2);
plot(V.t,sqrt(sum(V.Y.fSurface.^2,2)),'LineWidth',1.1); grid on;
ylabel('|d_F,m,surf| [N]'); xlabel('t [s]');
title('Thrust lag and actual omitted surface-force coupling');

figure;
plot(N.t,N.Y.w(:,1),'k--'); hold on;
plot(N.t,N.Y.what(:,1),'LineWidth',1.0);
plot(N.t,N.Y.w(:,2),'k--');
plot(N.t,N.Y.what(:,2),'LineWidth',1.0);
plot(N.t,N.Y.w(:,3),'k--');
plot(N.t,N.Y.what(:,3),'LineWidth',1.0);
grid on; xlabel('t [s]'); ylabel('wind [m/s]');
legend('w_N','hat w_N','w_E','hat w_E','w_D','hat w_D');
title('Noisy full-6DOF wind estimation');

figure;
plot3(V.Y.rd(:,1),V.Y.rd(:,2),-V.Y.rd(:,3),'k--','LineWidth',1.0); hold on;
plot3(V.Y.r(:,1),V.Y.r(:,2),-V.Y.r(:,3),'LineWidth',1.2);
grid on; axis equal;
xlabel('N [m]'); ylabel('E [m]'); zlabel('altitude [m]');
legend('r_d','r');
title('Complete nonlinear 6DOF ground-trajectory tracking');
end
