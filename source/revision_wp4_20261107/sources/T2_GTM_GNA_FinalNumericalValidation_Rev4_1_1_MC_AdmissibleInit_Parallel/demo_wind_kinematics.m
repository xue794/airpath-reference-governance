function demo_wind_kinematics()
%DEMO_WIND_KINEMATICS Demonstrate constant-wind Galilean equivalence.
P=t2_gna_params(); TR=t2_trim_level(P.src.V_mps,P.src.h_m,P);
w0=[8;-4;0]; U0=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
UW=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',w0);
phi0=deg2rad(2); theta0=TR.alpha+deg2rad(1); psi0=0; q0=euler321_to_quat(phi0,theta0,psi0);
va0=[TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha)]; C0=quat_to_dcm(q0);
x0=[0;0;-TR.h;va0;q0;0;0;0]; xw=x0; xw(4:6)=va0+C0.'*w0;
tg=linspace(0,25,501); opts=odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,X0]=ode45(@(t,x)t2_gna_6dof_quat(t,x,U0,P),tg,x0,opts);
[~,XW]=ode45(@(t,x)t2_gna_6dof_quat(t,x,UW,P),tg,xw,opts);
Va0=zeros(size(t)); VaW=zeros(size(t)); VgW=zeros(size(t)); posRes=zeros(size(t));
for k=1:numel(t)
    A0=t2_air_relative_state(t(k),X0(k,:).',U0,true); AW=t2_air_relative_state(t(k),XW(k,:).',UW,true);
    Va0(k)=A0.V; VaW(k)=AW.V; VgW(k)=norm(AW.vgB); posRes(k)=norm(XW(k,1:3).'-X0(k,1:3).'-w0*t(k));
end
figure; plot(t,Va0,'LineWidth',1.2); hold on; plot(t,VaW,'--','LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('airspeed V [m/s]'); legend('zero wind','constant wind');
figure; plot(t,VgW,'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('ground-speed body magnitude [m/s]');
figure; semilogy(t,max(posRes,eps),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('Galilean position residual [m]');
end
