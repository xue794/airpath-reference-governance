function O=viC_filter_error(Cc,omegaC,Cr,omegaR)
%VIC_FILTER_ERROR SO(3) command-filter error variables.

Q=Cc.'*Cr;
e=0.5*viC_vee(Q-Q.');
Psi=0.5*trace(eye(3)-Q);
xi=omegaR-Q.'*omegaC;

O.Q=Q;
O.e=e;
O.Psi=Psi;
O.xi=xi;
O.E=viC_E(Q);
O.theta=dcm_distance(Cc,Cr);
end
