function qdot = quat_kinematics(q,omega,kNorm)
%QUAT_KINEMATICS Quaternion kinematics for Cdot_IB=C_IB*hat(omega).
% Scalar-first q maps body to inertial. omega is body angular velocity.
% Optional kNorm >= 0 adds a norm-stabilization term that vanishes on S^3.
if nargin < 3, kNorm=0; end
q=q(:); omega=omega(:);
w=q(1); v=q(2:4);
qdot=0.5*[-v.'*omega; w*omega + cross(v,omega)];
if kNorm ~= 0
    qdot=qdot + kNorm*(1-q.'*q)*q;
end
end
