function xdot=vii_joint_rhs(t,x,Up,Cmd,R,A,F,Dp,P7)
%VII_JOINT_RHS VI-C command filter + physical T-2 moment allocation/servos.
%
% State:
%   [q_r(4); omega_r(3); q(4); omega(3); delta_act(3); T_act]
%
% Translational air-path states remain the verified upstream VI-A trajectory.
% Surface-induced translational forces are measured as d_F,m but are not
% fed into translational dynamics until the later full 6DOF integration.

qr=x(1:4);
omegaR=x(5:7);
q=x(8:11);
omega=x(12:14);
delta=x(15:17);
Tact=x(18); %#ok<NASGU>

Cc=viC_interp_command(t,Cmd);
Cr=quat_to_dcm(qr);
C=quat_to_dcm(q);

% VI-C command filter.
Qf=Cc.'*Cr;
eRf=0.5*viC_vee(Qf-Qf.');
qrDot=quat_kinematics(qr,omegaR,F.quatNormGain);
omegaRDot=-F.kRf*eRf-F.KomegaF*omegaR;

% Prescribed upstream path direction and speed.
V=interp1(Up.t,Up.y(:,1),t,'pchip');
gamma=interp1(Up.t,Up.y(:,2),t,'pchip');
chi=interp1(Up.t,Up.y(:,3),t,'pchip');
Fp=interp1(Up.t,Up.Fparallel,t,'pchip');

Pth=viA_airpath_basis(gamma,chi);
Air=viB_airdata_from_geometry(C,Pth.t);

% VI-D virtual torque with the same physical D_omega,a.
Mom=vii_moment_components(V,Air.alpha,Air.beta,R);
Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot,C,omega,Mom.D,R,Dp);

% Constrained allocation.
Alloc=vii_allocate(Ctl.tauC,V,Air.alpha,Air.beta,omega,R,P7);

% Body-axis thrust correction uses the ACTUAL alpha,beta geometry.
Tcmd=Fp/Air.kappa;

Act=vii_actuator_derivative(delta,x(18),Alloc.deltaCmd,Tcmd,P7);

% Exact physical GNA moment model.
tauA=Mom.tau0+Mom.K*omega+Mom.B*delta;
omegaDot=R.P.I\(tauA-cross(omega,R.P.I*omega));
qDot=quat_kinematics(q,omega,F.quatNormGain);

xdot=[qrDot;omegaRDot;qDot;omegaDot;Act.deltaDot;Act.TDot];
end
