function xdot = t2_gna_6dof_euler(t,x,U,P)
%T2_GNA_6DOF_EULER  Open-loop nonlinear 6DOF plant in NED/body FRD axes.
%
% State x = [N E D u_g v_g w_g phi theta psi wx wy wz]'.
% Controls U.delta=[da de dr]' [rad], U.T [N], U.wI=[wN wE wD]' [m/s].
% vgB is ground-relative velocity expressed in body axes; aerodynamics use
% vaB = vgB - C_IB' wI.

rI=x(1:3); %#ok<NASGU>
vgB=x(4:6); phi=x(7); theta=x(8); psi=x(9); omega=x(10:12);
CIB=euler321_dcm(phi,theta,psi);
[wI,~]=t2_wind_eval(t,U);
h = -x(3);  % NED D positive down, altitude positive up
[fA,tauA] = t2_gna_aero(vgB,omega,U.delta,CIB,wI,h,P);

rDot = CIB*vgB;
zeta = CIB.'*[0;0;1];
fT = [U.T;0;0];
vgDot = -cross(omega,vgB) + P.g*zeta + (fA+fT)/P.m;
omegaDot = P.I \ (tauA - cross(omega,P.I*omega));

p=omega(1); q=omega(2); r=omega(3);
if abs(cos(theta)) < 1e-6
    error('t2_gna_6dof_euler:EulerSingularity','theta too close to +/-90 deg.');
end
phiDot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
thetaDot = q*cos(phi)-r*sin(phi);
psiDot = (q*sin(phi)+r*cos(phi))/cos(theta);
xdot=[rDot;vgDot;phiDot;thetaDot;psiDot;omegaDot];
end
