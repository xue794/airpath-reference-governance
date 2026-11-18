function demo_section_viA(mode)
%DEMO_SECTION_VIA Diagnostic plots for Section VI-A.
if nargin<1, mode='varying'; end

R=rg_params(); A=viA_params(); E=wind_estimator_params();
switch lower(mode)
    case 'still', T=A.Tstill;
    case 'constant', T=A.Tconstant;
    case 'varying', T=A.Tvarying;
    case 'aggressive', T=A.Taggressive;
    case 'wrap', T=A.Twrap;
    otherwise, error('Unknown mode');
end

S=simulate_viA_case(mode,T,R,A,E);

figure;
subplot(3,1,1);
plot(S.t,S.xf(:,1),'k--','LineWidth',1.1); hold on;
plot(S.t,S.y(:,1),'LineWidth',1.2); grid on;
ylabel('V [m/s]'); legend('V_f','V');
subplot(3,1,2);
plot(S.t,rad2deg(S.xf(:,2)),'k--','LineWidth',1.1); hold on;
plot(S.t,rad2deg(S.y(:,2)),'LineWidth',1.2); grid on;
ylabel('\gamma [deg]'); legend('\gamma_f','\gamma');
subplot(3,1,3);
plot(S.t,rad2deg(S.xf(:,3)),'k--','LineWidth',1.1); hold on;
plot(S.t,rad2deg(S.y(:,3)),'LineWidth',1.2); grid on;
ylabel('\chi [deg]'); xlabel('t [s]'); legend('\chi_f','\chi');

figure;
subplot(3,1,1);
plot(S.t,S.e(:,1),'LineWidth',1.2); grid on; ylabel('e_V [m/s]');
subplot(3,1,2);
plot(S.t,rad2deg(S.e(:,2)),'LineWidth',1.2); grid on; ylabel('e_\gamma [deg]');
subplot(3,1,3);
plot(S.t,rad2deg(S.e(:,3)),'LineWidth',1.2); grid on;
ylabel('e_\chi [deg]'); xlabel('t [s]');

figure;
subplot(3,1,1);
plot(S.t,rad2deg(S.alpha),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(S.alphaFF),'k--','LineWidth',1.0);
plot(S.t,rad2deg(R.alphaLo)*ones(size(S.t)),'k:');
plot(S.t,rad2deg(R.alphaHi)*ones(size(S.t)),'k:');
grid on; ylabel('\alpha [deg]'); legend('\alpha_c','\alpha_{ff}');
subplot(3,1,2);
plot(S.t,S.Tstar,'LineWidth',1.2); hold on;
plot(S.t,S.Tff,'k--','LineWidth',1.0);
plot(S.t,R.Tmin*ones(size(S.t)),'k:');
plot(S.t,R.Tmax*ones(size(S.t)),'k:');
grid on; ylabel('T [N]'); legend('T^*','T_{ff}');
subplot(3,1,3);
plot(S.t,rad2deg(S.mu),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(S.muFF),'k--','LineWidth',1.0);
plot(S.t,rad2deg(R.muMax)*ones(size(S.t)),'k:');
plot(S.t,-rad2deg(R.muMax)*ones(size(S.t)),'k:');
grid on; ylabel('\mu [deg]'); xlabel('t [s]'); legend('\mu_c','\mu_{ff}');

if strcmpi(mode,'varying')
    figure;
    plot(S.t,S.wdot,'LineWidth',1.2); grid on;
    xlabel('t [s]'); ylabel('wind acceleration [m/s^2]');
    legend('dw_N/dt','dw_E/dt','dw_D/dt');
    title('Time-varying wind disturbance entering ideal air-path errors');
end
end
