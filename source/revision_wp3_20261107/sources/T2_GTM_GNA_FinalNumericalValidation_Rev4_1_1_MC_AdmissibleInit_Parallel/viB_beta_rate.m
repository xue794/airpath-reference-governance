function betaDot=viB_beta_rate(C,t,omega,tdot)
%VIB_BETA_RATE Exact sideslip derivative identity.
%
% bdot = b x omega + C'*tdot
% cos(beta)*betaDot = e2'*bdot

A=viB_airdata_from_geometry(C,t);
bdot=cross(A.b,omega)+C.'*tdot;
betaDot=bdot(2)/cos(A.beta);
end
