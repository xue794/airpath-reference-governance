function xdot=viD_joint_rhs(t,x,Up,Cmd,distMode,R,A,F,Dp)
%VID_JOINT_RHS Dynamic command filter + ideal rotational PHS plant.
%
% State:
%   x=[q_r(4); omega_r(3); q(4); omega(3)]
%
% Upstream Section V/VI-A is one-way into VI-C/VI-D and is therefore
% precomputed and interpolated.

qr=x(1:4);
omegaR=x(5:7);
q=x(8:11);
omega=x(12:14);

Cc=viC_interp_command(t,Cmd);
Cr=quat_to_dcm(qr);
C=quat_to_dcm(q);

% Command filter.
Qf=Cc.'*Cr;
eRf=0.5*viC_vee(Qf-Qf.');
qrDot=quat_kinematics(qr,omegaR,F.quatNormGain);
omegaRDot=-F.kRf*eRf-F.KomegaF*omegaR;

% Actual air-path direction, used only to evaluate D_omega,a(V,alpha).
V=interp1(Up.t,Up.y(:,1),t,'pchip');
gamma=interp1(Up.t,Up.y(:,2),t,'pchip');
chi=interp1(Up.t,Up.y(:,3),t,'pchip');
Bpath=viA_airpath_basis(gamma,chi);
Air=viB_airdata_from_geometry(C,Bpath.t);

Mom=viD_rate_moment_decomp(V,Air.alpha,R);
Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot,C,omega,Mom.D,R,Dp);

dTau=viD_torque_disturbance(t,distMode,Dp);

omegaDot=R.P.I\(-cross(omega,R.P.I*omega)-Mom.D*omega+Ctl.tauC+dTau);
qDot=quat_kinematics(q,omega,F.quatNormGain);

xdot=[qrDot;omegaRDot;qDot;omegaDot];
end
