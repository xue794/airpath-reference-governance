function K=viC_command_kinematics(t,z,mode,R,A,Ew)
%VIC_COMMAND_KINEMATICS Analytic command-attitude angular velocity.
%
% This function uses the current governor acceleration u*=ddot(c_f), but it
% never differentiates u*.  Hence omega_c requires no hidden third
% derivative of the governed air-path reference.

xf=z(1:6);
y=z(7:9);

Sig=viA_case_signal(t,mode,Ew);
[vc,~]=rg_raw_test_command(t,Sig.refMode,R);
[cc,~]=rg_raw_command_chart(vc,R);
u0=rg_nominal_control(xf,cc,R);
Gov=rg_scalar_governor(xf,u0,R);
xfDot=[xf(4:6);Gov.u];

[yDot,C]=viA_ideal_path_rhs(y,xf,Sig.wdot,R,A);
eDot=yDot-xf(4:6);

V=y(1); gamma=y(2);
nuG=xf(5); nuC=xf(6);

qGamma=nuG-A.kGamma*C.e(2);
qChi=nuC-A.kChi*C.e(3);

Xdot=R.P.m*(Gov.u(1)-A.kV*eDot(1)+R.g*cos(gamma)*yDot(2));

aGammaDot= ...
    yDot(1)*qGamma + ...
    V*(Gov.u(2)-A.kGamma*eDot(2)) - ...
    R.g*sin(gamma)*yDot(2);

aChiDot= ...
    yDot(1)*cos(gamma)*qChi - ...
    V*sin(gamma)*yDot(2)*qChi + ...
    V*cos(gamma)*(Gov.u(3)-A.kChi*eDot(3));

rr=hypot(C.aGamma,C.aChi);
Ndot=R.P.m*(C.aGamma*aGammaDot+C.aChi*aChiDot)/rr;

muDot=(C.aGamma*aChiDot-C.aChi*aGammaDot) / ...
      (C.aGamma^2+C.aChi^2);

M=rg_N_map(C.alpha,V,C.X,R);
alphaDot=(Ndot-M.V*yDot(1)-M.X*Xdot)/M.alpha;

a=C.alpha;
mu=C.mu;

e1=[1;0;0]; e2=[0;1;0]; e3=[0;0;1];

omegaC= ...
    yDot(3)*(viB_Ry(a).'*viB_Rx(mu).'*viB_Ry(gamma).'*e3) + ...
    yDot(2)*(viB_Ry(a).'*viB_Rx(mu).'*e2) + ...
    muDot*(viB_Ry(a).'*e1) + ...
    alphaDot*e2;

K.C=C;
K.Cc=viB_command_attitude(gamma,y(3),mu,a);
K.omegaC=omegaC;
K.alphaDot=alphaDot;
K.muDot=muDot;
K.Xdot=Xdot;
K.Ndot=Ndot;
K.aGammaDot=aGammaDot;
K.aChiDot=aChiDot;
K.upDot=[xfDot;yDot];
K.u=Gov.u;
K.lambda=Gov.lambda;
K.yDot=yDot;
K.eDot=eDot;
end
