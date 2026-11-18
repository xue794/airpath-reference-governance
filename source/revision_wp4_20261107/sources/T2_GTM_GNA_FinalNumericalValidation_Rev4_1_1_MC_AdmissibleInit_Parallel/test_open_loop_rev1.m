function Results = test_open_loop_rev1()
%TEST_OPEN_LOOP_REV1 P1--P6 validation for quaternion/SO(3)+wind Rev.1.
% MATLAB R2018a compatible; base MATLAB only.
R0=test_open_loop(); P=R0.P; TR=R0.trim;
fprintf('\nT-2 GTM--GNA Open-Loop Rev.1 validation\n');
fprintf('======================================\n');

%% P4: Euler vs quaternion equivalence in still air
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
phi0=deg2rad(2); theta0=TR.alpha+deg2rad(1); psi0=0;
vg0=[TR.V*cos(TR.alpha);0;TR.V*sin(TR.alpha)];
xE0=[0;0;-TR.h;vg0;phi0;theta0;psi0;0;0;0];
q0=euler321_to_quat(phi0,theta0,psi0);
xQ0=[0;0;-TR.h;vg0;q0;0;0;0];
tg=linspace(0,40,801); opts=odeset('RelTol',1e-10,'AbsTol',1e-11);
[tE,xE]=ode45(@(t,x)t2_gna_6dof_euler(t,x,U,P),tg,xE0,opts);
[tQ,xQ]=ode45(@(t,x)t2_gna_6dof_quat(t,x,U,P),tg,xQ0,opts); %#ok<ASGLU>
errR=0; errV=0; errW=0; errC=0; qNormErr=0;
for k=1:numel(tE)
    CE=euler321_dcm(xE(k,7),xE(k,8),xE(k,9));
    CQ=quat_to_dcm(xQ(k,7:10).');
    errR=max(errR,norm(xE(k,1:3).'-xQ(k,1:3).'));
    errV=max(errV,norm(xE(k,4:6).'-xQ(k,4:6).'));
    errW=max(errW,norm(xE(k,10:12).'-xQ(k,11:13).'));
    errC=max(errC,dcm_distance(CE,CQ));
    qNormErr=max(qNormErr,abs(norm(xQ(k,7:10))-1));
end
fprintf('P4 Euler/quaternion equivalence (40 s):\n');
fprintf('   max position error = %.3e m\n',errR);
fprintf('   max velocity error = %.3e m/s\n',errV);
fprintf('   max rate error = %.3e rad/s\n',errW);
fprintf('   max attitude geodesic error = %.3e rad\n',errC);
fprintf('   max quaternion norm error = %.3e\n',qNormErr);
assert(errR<1e-5 && errV<1e-7 && errW<1e-7 && errC<1e-7 && qNormErr<1e-8,'P4 failed');

%% P5: constant horizontal wind Galilean-equivalence test
w0=[8;-4;0];
U0=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
UW=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',w0);
% Use the same perturbed attitude and same initial AIR-relative body velocity.
C0=quat_to_dcm(q0);
xQ0a=xQ0;
xQ0w=xQ0; xQ0w(4:6)=vg0+C0.'*w0;
tg2=linspace(0,25,501);
[t0,x0]=ode45(@(t,x)t2_gna_6dof_quat(t,x,U0,P),tg2,xQ0a,opts);
[tw,xw]=ode45(@(t,x)t2_gna_6dof_quat(t,x,UW,P),tg2,xQ0w,opts); %#ok<ASGLU>
errAir=0; errAtt=0; errRate=0; errGal=0;
for k=1:numel(t0)
    A0=t2_air_relative_state(t0(k),x0(k,:).',U0,true);
    AW=t2_air_relative_state(t0(k),xw(k,:).',UW,true);
    errAir=max(errAir,norm(A0.vaB-AW.vaB));
    errAtt=max(errAtt,dcm_distance(A0.CIB,AW.CIB));
    errRate=max(errRate,norm(A0.omega-AW.omega));
    dpos=xw(k,1:3).'-x0(k,1:3).'-w0*t0(k);
    errGal=max(errGal,norm(dpos));
end
fprintf('P5 constant-wind Galilean equivalence, w=[8 -4 0] m/s:\n');
fprintf('   max air-relative velocity mismatch = %.3e m/s\n',errAir);
fprintf('   max attitude mismatch = %.3e rad\n',errAtt);
fprintf('   max angular-rate mismatch = %.3e rad/s\n',errRate);
fprintf('   max position Galilean residual = %.3e m\n',errGal);
assert(errAir<1e-7 && errAtt<1e-7 && errRate<1e-7 && errGal<1e-5,'P5 failed');

%% P6: time-varying wind derivative identity
Utv=struct('delta',[deg2rad(1.5);TR.deltae+deg2rad(0.4);deg2rad(-2)], ...
    'T',TR.T+2.0,'wind_fun',@wind_profile);
tstar=3.7;
qstar=euler321_to_quat(deg2rad(8),deg2rad(4),deg2rad(20));
Cstar=quat_to_dcm(qstar);
[wstar,wdstar]=t2_wind_eval(tstar,Utv);
vastar=[35;1.2;3.1];
vgstar=vastar+Cstar.'*wstar;
xstar=[120;-30;-350;vgstar;qstar;deg2rad([6;-3;4])];
xd=t2_gna_6dof_quat(tstar,xstar,Utv,P);
Astar=t2_air_relative_state(tstar,xstar,Utv,true);
[fA,~]=t2_gna_aero(Astar.vgB,Astar.omega,Utv.delta,Astar.CIB,wstar,-xstar(3),P);
zeta=Astar.CIB.'*[0;0;1]; fT=[Utv.T;0;0];
vaDotFromGround=xd(4:6)+cross(Astar.omega,Astar.CIB.'*wstar)-Astar.CIB.'*wdstar;
vaDotDirect=-cross(Astar.omega,Astar.vaB)+P.g*zeta+(fA+fT)/P.m-Astar.CIB.'*wdstar;
errWindDot=norm(vaDotFromGround-vaDotDirect);
fprintf('P6 time-varying-wind air-relative derivative identity:\n');
fprintf('   ||vaDot(ground identity)-vaDot(direct)|| = %.3e m/s^2\n',errWindDot);
assert(errWindDot<1e-11,'P6 failed');

Results=struct('Rev0',R0,'P4',struct('pos',errR,'vel',errV,'rate',errW,'att',errC,'qnorm',qNormErr), ...
    'P5',struct('air',errAir,'att',errAtt,'rate',errRate,'galilean',errGal), ...
    'P6',struct('identity',errWindDot));

    function z=wind_profile(t)
        w=[4+2*sin(0.4*t); -3+1.5*cos(0.7*t); 0.5*sin(0.3*t)];
        wd=[0.8*cos(0.4*t); -1.05*sin(0.7*t); 0.15*cos(0.3*t)];
        z=[w;wd];
    end
end
