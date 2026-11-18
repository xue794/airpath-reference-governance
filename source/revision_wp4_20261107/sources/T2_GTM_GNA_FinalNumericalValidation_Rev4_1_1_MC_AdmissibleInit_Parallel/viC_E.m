function E=viC_E(Q)
%VIC_E Exact attitude-error kinematics matrix for this paper's Q convention.
%
% Q = C_c' C_r
% e = 0.5*(Q-Q')^vee
% Qdot = Q*hat(xi)
%
% Therefore
%   edot = E(Q)*xi
% with
%   E(Q)=0.5*(trace(Q)I-Q').
%
% NOTE: the transpose is essential for the present Q convention.
E=0.5*(trace(Q)*eye(3)-Q.');
end
