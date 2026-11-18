function demo_open_loop()
%DEMO_OPEN_LOOP Small-disturbance open-loop simulation around level trim.
P=t2_gna_params(); TR=t2_trim_level(P.src.V_mps,P.src.h_m,P);
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
% State [N E D u_g v_g w_g phi theta psi wx wy wz]
x0=[0;0;-TR.h; TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha); ...
    deg2rad(2);TR.alpha+deg2rad(1);0; 0;0;0];
opts=odeset('RelTol',1e-9,'AbsTol',1e-10);
[t,x]=ode45(@(t,x)t2_gna_6dof_euler(t,x,U,P),[0 40],x0,opts);
V=zeros(size(t)); alpha=zeros(size(t));
for k=1:numel(t)
    A=t2_air_relative_state(t(k),x(k,:).',U,false); V(k)=A.V; alpha(k)=A.alpha;
end
figure; plot(t,V,'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('Airspeed V [m/s]');
figure; plot(t,rad2deg(alpha),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('Angle of attack \alpha [deg]');
figure; plot(t,rad2deg(x(:,7:9)),'LineWidth',1.2); grid on; xlabel('t [s]'); ylabel('Euler angles [deg]'); legend('\phi','\theta','\psi');
end
