function xdot=viD_constant_rhs(~,x,Cr,V,alpha,R,F,Dp)
%VID_CONSTANT_RHS Fixed-reference ideal rotational PHS test.
%
% State x=[q(4);omega(3)], with omega_r=omegaDot_r=0.

q=x(1:4);
omega=x(5:7);
C=quat_to_dcm(q);

Mom=viD_rate_moment_decomp(V,alpha,R);
Ctl=viD_virtual_torque(Cr,zeros(3,1),zeros(3,1),C,omega,Mom.D,R,Dp);

omegaDot=R.P.I\(-cross(omega,R.P.I*omega)-Mom.D*omega+Ctl.tauC);
qDot=quat_kinematics(q,omega,F.quatNormGain);

xdot=[qDot;omegaDot];
end
