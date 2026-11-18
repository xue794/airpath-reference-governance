function xdot = t2_reduced9_dynamics(x,U,h,P)
%T2_REDUCED9_DYNAMICS Validation-only dynamics at fixed density altitude.
% x=[u v w wx wy wz phi theta psi]'. Position is omitted.
vgB=x(1:3); omega=x(4:6); phi=x(7); theta=x(8); psi=x(9);
CIB=euler321_dcm(phi,theta,psi);
[wI,~]=t2_wind_eval(0,U);
[fA,tauA]=t2_gna_aero(vgB,omega,U.delta,CIB,wI,h,P);
zeta=CIB.'*[0;0;1];
vgDot=-cross(omega,vgB)+P.g*zeta+(fA+[U.T;0;0])/P.m;
omegaDot=P.I\(tauA-cross(omega,P.I*omega));
p=omega(1); q=omega(2); r=omega(3);
phiDot=p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
thetaDot=q*cos(phi)-r*sin(phi);
psiDot=(q*sin(phi)+r*cos(phi))/cos(theta);
xdot=[vgDot;omegaDot;phiDot;thetaDot;psiDot];
end
