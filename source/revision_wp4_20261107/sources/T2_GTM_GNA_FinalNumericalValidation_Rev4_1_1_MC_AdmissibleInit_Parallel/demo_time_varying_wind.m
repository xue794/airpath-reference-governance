function demo_time_varying_wind()
%DEMO_TIME_VARYING_WIND Open-loop quaternion plant under smooth uniform wind.
P=t2_gna_params(); TR=t2_trim_level(P.src.V_mps,P.src.h_m,P);
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wind_fun',@t2_demo_wind_profile);
q0=euler321_to_quat(0,TR.alpha,0);
C0=quat_to_dcm(q0); w0=t2_demo_wind_profile(0); w0=w0(1:3);
va0=[TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha)];
vg0=va0+C0.'*w0; % initialize at the same air-relative trim state
x0=[0;0;-TR.h;vg0;q0;0;0;0];
tg=linspace(0,40,801); opts=odeset('RelTol',1e-9,'AbsTol',1e-10);
[t,x]=ode45(@(t,x)t2_gna_6dof_quat(t,x,U,P),tg,x0,opts);
Va=zeros(size(t)); Vg=zeros(size(t)); alpha=zeros(size(t)); W=zeros(numel(t),3);
for k=1:numel(t)
    A=t2_air_relative_state(t(k),x(k,:).',U,true);
    Va(k)=A.V; Vg(k)=norm(A.vgB); alpha(k)=A.alpha; W(k,:)=A.wI.';
end
figure; plot(t,Va,'LineWidth',1.2); hold on; plot(t,Vg,'--','LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('speed [m/s]'); legend('airspeed V','ground-speed magnitude');
figure; plot(t,rad2deg(alpha),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('angle of attack \alpha [deg]');
figure; plot(t,W,'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('wind components [m/s]'); legend('w_N','w_E','w_D');
end
