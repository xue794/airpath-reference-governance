function O=viD_virtual_torque(Cr,omegaR,omegaRDot,C,omega,Domega,R,Dp)
%VID_VIRTUAL_TORQUE PHS-based rotational virtual torque, Eq. (98) convention.
%
% tau_c =
%   omega x I omega
%   + D_omega,a A_R omega_r
%   - k_R e_R
%   - K_omega,c e_omega
%   - I hat(omega) A_R omega_r
%   + I A_R omegaDot_r.
%
% The ideal VI-D plant is
%   I omegaDot + omega x I omega = -D_omega,a omega + tau_c + d_tau.
%
% Thus
%   I eOmegaDot = -k_R e_R -(K_omega,c+D_omega,a)eOmega + d_tau.

Err=viD_error(Cr,omegaR,C,omega);
I=R.P.I;

tauC=cross(omega,I*omega) + ...
     Domega*(Err.A*omegaR) - ...
     Dp.kR*Err.eR - ...
     Dp.KomegaC*Err.eOmega - ...
     I*viC_hat(omega)*(Err.A*omegaR) + ...
     I*(Err.A*omegaRDot);

O.tauC=tauC;
O.Err=Err;
O.Kd=Dp.KomegaC+Domega;
end
