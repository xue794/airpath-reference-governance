function xdot = t2_gna_6dof_quat(t,x,U,P)
%T2_GNA_6DOF_QUAT Open-loop nonlinear 6DOF plant, quaternion/SO(3) form.
%
% State x = [N E D u_g v_g w_g q0 q1 q2 q3 wx wy wz]'.
% q is scalar-first and maps body B -> NED inertial I.
% Ground-relative velocity vgB is expressed in body axes. Aerodynamics use
% vaB = vgB - C_IB' wI. Spatially uniform wind is evaluated in inertial axes.

vgB=x(4:6); qIB=x(7:10); omega=x(11:13);
CIB=quat_to_dcm(qIB);
[wI,~]=t2_wind_eval(t,U);
h=-x(3);
[fA,tauA]=t2_gna_aero(vgB,omega,U.delta,CIB,wI,h,P);

rDot=CIB*vgB;
zeta=CIB.'*[0;0;1];
fT=[U.T;0;0];
vgDot=-cross(omega,vgB)+P.g*zeta+(fA+fT)/P.m;
omegaDot=P.I\(tauA-cross(omega,P.I*omega));

if isfield(P,'quat_norm_gain')
    kNorm=P.quat_norm_gain;
else
    kNorm=0;
end
qDot=quat_kinematics(qIB,omega,kNorm);
xdot=[rDot;vgDot;qDot;omegaDot];
end
