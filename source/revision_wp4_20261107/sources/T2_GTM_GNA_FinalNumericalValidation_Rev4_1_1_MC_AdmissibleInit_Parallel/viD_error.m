function E=viD_error(Cr,omegaR,C,omega)
%VID_ERROR Rotational tracking errors.
%
% Q_R = C_r' C
% A_R = C' C_r = Q_R'
% e_R = 0.5(Q_R-Q_R')^vee
% e_omega = omega-A_R*omega_r

Q=Cr.'*C;
Arel=Q.';
eR=0.5*viC_vee(Q-Q.');
eOmega=omega-Arel*omegaR;
Psi=0.5*trace(eye(3)-Q);

E.Q=Q;
E.A=Arel;
E.eR=eR;
E.eOmega=eOmega;
E.Psi=Psi;
E.theta=dcm_distance(Cr,C);
E.E=viC_E(Q);
end
