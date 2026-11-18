function demo_quat_equivalence()
%DEMO_QUAT_EQUIVALENCE Compare Euler and quaternion plants in still air.
P=t2_gna_params(); TR=t2_trim_level(P.src.V_mps,P.src.h_m,P);
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
phi0=deg2rad(2); theta0=TR.alpha+deg2rad(1); psi0=0;
vg0=[TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha)];
xE0=[0;0;-TR.h;vg0;phi0;theta0;psi0;0;0;0];
xQ0=[0;0;-TR.h;vg0;euler321_to_quat(phi0,theta0,psi0);0;0;0];
tg=linspace(0,40,801); opts=odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,xE]=ode45(@(t,x)t2_gna_6dof_euler(t,x,U,P),tg,xE0,opts);
[~,xQ]=ode45(@(t,x)t2_gna_6dof_quat(t,x,U,P),tg,xQ0,opts);
attErr=zeros(size(t)); velErr=zeros(size(t)); qNorm=zeros(size(t));
for k=1:numel(t)
    CE=euler321_dcm(xE(k,7),xE(k,8),xE(k,9)); CQ=quat_to_dcm(xQ(k,7:10).');
    attErr(k)=dcm_distance(CE,CQ); velErr(k)=norm(xE(k,4:6)-xQ(k,4:6)); qNorm(k)=norm(xQ(k,7:10));
end
figure; semilogy(t,max(attErr,eps),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('SO(3) attitude mismatch [rad]');
figure; semilogy(t,max(velErr,eps),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('velocity mismatch [m/s]');
figure; plot(t,qNorm,'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('quaternion norm');
end
