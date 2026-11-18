function demo_reference_governor(mode)
%DEMO_REFERENCE_GOVERNOR Section V reference-governor diagnostic figures.
% MATLAB R2018a compatible (no yline/tiledlayout).
if nargin<1, mode='aggressive'; end
R=rg_params();
if strcmpi(mode,'small'), T=R.Tsmall;
elseif strcmpi(mode,'wrap'), T=R.Twrap;
else, T=R.Taggressive;
end
S=simulate_reference_governor(mode,T,R);

figure;
subplot(3,1,1);
plot(S.t,S.cc(:,1),'k--','LineWidth',1.1); hold on;
plot(S.t,S.x(:,1),'LineWidth',1.2); grid on;
ylabel('V [m/s]'); legend('V_c','V_f');
subplot(3,1,2);
plot(S.t,rad2deg(S.cc(:,2)),'k--','LineWidth',1.1); hold on;
plot(S.t,rad2deg(S.x(:,2)),'LineWidth',1.2); grid on;
ylabel('\gamma [deg]'); legend('\gamma_c','\gamma_f');
subplot(3,1,3);
plot(S.t,rad2deg(unwrap(S.cc(:,3))),'k--','LineWidth',1.1); hold on;
plot(S.t,rad2deg(S.x(:,3)),'LineWidth',1.2); grid on;
ylabel('\chi [deg]'); xlabel('t [s]'); legend('\chi_c','\chi_f');

figure;
plot(S.t,S.lambda,'LineWidth',1.2); grid on; ylim([0 1.05]);
xlabel('t [s]'); ylabel('\lambda^*'); title('Certified scalar governor');

figure;
plot(S.t,S.minphi,'LineWidth',1.2); hold on;
plot(S.t,zeros(size(S.t)),'k--'); grid on;
xlabel('t [s]'); ylabel('min_j \phi_j'); title('Minimum CBF/HOCBF feasibility margin');

figure;
subplot(3,1,1);
plot(S.t,rad2deg(S.alpha),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(R.alphaLo)*ones(size(S.t)),'k--');
plot(S.t,rad2deg(R.alphaHi)*ones(size(S.t)),'k--'); grid on;
ylabel('\alpha_{ff} [deg]');
subplot(3,1,2);
plot(S.t,S.Tff,'LineWidth',1.2); hold on;
plot(S.t,(R.Tmin+R.DeltaT)*ones(size(S.t)),'k--');
plot(S.t,(R.Tmax-R.DeltaT)*ones(size(S.t)),'k--'); grid on;
ylabel('T_{ff} [N]');
subplot(3,1,3);
plot(S.t,rad2deg(S.mu),'LineWidth',1.2); hold on;
plot(S.t,rad2deg(R.muMax)*ones(size(S.t)),'k--');
plot(S.t,-rad2deg(R.muMax)*ones(size(S.t)),'k--'); grid on;
ylabel('\mu_{ff} [deg]'); xlabel('t [s]');

figure;
plot(S.t,S.hB,'LineWidth',1.2); hold on;
plot(S.t,zeros(size(S.t)),'k--'); grid on;
xlabel('t [s]'); ylabel('h_B'); title('Backup-set barrier');
end
