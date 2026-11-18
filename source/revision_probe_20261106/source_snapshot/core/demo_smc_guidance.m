function demo_smc_guidance(mode)
%DEMO_SMC_GUIDANCE Plot isolated guidance behavior.
% Example:
%   demo_smc_guidance('varying')
%   demo_smc_guidance('bound')
%   demo_smc_guidance('noisy')

if nargin<1, mode='varying'; end
E=wind_estimator_params(); G=smc_guidance_params();

if strcmpi(mode,'still') || strcmpi(mode,'constant')
    T=60;
elseif strcmpi(mode,'mismatch')
    T=80;
else
    T=140;
end

R=simulate_smc_guidance_case(mode,T,G,E);

figure;
plot(R.t,R.er,'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('e_r [m]');
legend('N','E','D'); title(['Position tracking error, mode = ' mode]);

figure;
plot(R.t,R.s,'LineWidth',1.2); hold on; grid on;
plot(R.t, G.eps(1)*ones(size(R.t)),'k--');
plot(R.t,-G.eps(1)*ones(size(R.t)),'k--');
xlabel('t [s]'); ylabel('s [m]');
legend('s_N','s_E','s_D','+eps_N','-eps_N');
title('Sliding variable');

figure;
plot(R.t,R.w(:,1),'--','LineWidth',1.2); hold on;
plot(R.t,R.what(:,1),'LineWidth',1.2);
plot(R.t,R.w(:,2),'--','LineWidth',1.2);
plot(R.t,R.what(:,2),'LineWidth',1.2);
plot(R.t,R.w(:,3),'--','LineWidth',1.2);
plot(R.t,R.what(:,3),'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('Wind [m/s]');
legend('w_N','hat w_N','w_E','hat w_E','w_D','hat w_D');
title('True and estimated wind');

figure;
plot3(R.rd(:,2),R.rd(:,1),-R.rd(:,3),'k--','LineWidth',1.2); hold on;
plot3(R.r(:,2),R.r(:,1),-R.r(:,3),'LineWidth',1.2); grid on; axis equal;
xlabel('E [m]'); ylabel('N [m]'); zlabel('Altitude [m]');
legend('r_d','r'); title('Ground-trajectory tracking');

if any(strcmpi(mode,{'mismatch','bound','noisy'}))
    figure;
    plot(R.t,R.eg,'LineWidth',1.2); grid on;
    xlabel('t [s]'); ylabel('e_g [m/s]');
    legend('N','E','D'); title('Synthetic governor mismatch');
end
end
