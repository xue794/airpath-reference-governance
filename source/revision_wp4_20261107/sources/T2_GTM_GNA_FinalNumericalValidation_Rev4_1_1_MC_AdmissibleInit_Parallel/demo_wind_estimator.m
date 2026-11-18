function demo_wind_estimator(addNoise)
%DEMO_WIND_ESTIMATOR Time-varying wind estimation on the verified quaternion plant.
if nargin<1, addNoise=true; end
P=t2_gna_params(); TR=t2_trim_level(P.src.V_mps,P.src.h_m,P); E=wind_estimator_params();
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wind_fun',@windfun);
[w0,~]=t2_wind_eval(0,U); q0=euler321_to_quat(0,TR.alpha,0); C0=quat_to_dcm(q0);
va0=[TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha)]; vg0=va0+C0.'*w0;
x0=[0;0;-TR.h;vg0;q0;0;0;0]; z0=[x0;E.what0];
tg=linspace(0,60,2401); opts=odeset('RelTol',1e-9,'AbsTol',1e-10);
[t,z]=ode45(@(t,z)wind_estimator_augmented_rhs(t,z,U,P,E,addNoise),tg,z0,opts);
W=zeros(numel(t),3); WH=z(:,14:16); Va=zeros(numel(t),1); Vg=zeros(numel(t),1);
for k=1:numel(t)
    [wk,~]=t2_wind_eval(t(k),U); W(k,:)=wk.';
    C=quat_to_dcm(z(k,7:10).'); vgB=z(k,4:6).'; vaB=vgB-C.'*wk;
    Va(k)=norm(vaB); Vg(k)=norm(vgB);
end
err=W-WH;
figure; plot(t,W,'--','LineWidth',1.2); hold on; plot(t,WH,'LineWidth',1.1); grid on;
xlabel('t [s]'); ylabel('Wind velocity [m/s]');
legend('w_N','w_E','w_D','hat w_N','hat w_E','hat w_D','Location','best');
title(sprintf('Wind estimator, measurement noise = %d',logical(addNoise)));
figure; plot(t,err,'LineWidth',1.1); grid on; xlabel('t [s]'); ylabel('w-hat w [m/s]');
legend('N','E','D','Location','best'); title('Wind estimation error');
figure; plot(t,Va,'LineWidth',1.1); hold on; plot(t,Vg,'LineWidth',1.1); grid on;
xlabel('t [s]'); ylabel('Speed [m/s]'); legend('Airspeed V_a','Ground speed V_g','Location','best');
end

function z=windfun(t)
[w,wd]=wind_estimator_test_wind(t,'varying'); z=[w;wd];
end
