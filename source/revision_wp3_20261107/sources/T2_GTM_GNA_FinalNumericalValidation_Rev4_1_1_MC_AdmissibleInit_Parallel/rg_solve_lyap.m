function P = rg_solve_lyap(A,Q)
%RG_SOLVE_LYAP Solve A'*P + P*A = -Q using base MATLAB only.
% Avoids dependency on Control System Toolbox lyap().
n=size(A,1);
K = kron(eye(n),A') + kron(A',eye(n));
p = -K\Q(:);
P = reshape(p,n,n);
P = 0.5*(P+P');
end
