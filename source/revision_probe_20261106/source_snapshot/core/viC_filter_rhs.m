function xdot=viC_filter_rhs(~,x,Cc,F)
%VIC_FILTER_RHS Command-filter dynamics for a supplied command DCM.
%
% State x=[q_r(4);omega_r(3)].

qr=x(1:4);
omegaR=x(5:7);
Cr=quat_to_dcm(qr);

Err=viC_filter_error(Cc,zeros(3,1),Cr,omegaR);

qrDot=quat_kinematics(qr,omegaR,F.quatNormGain);
omegaRDot=-F.kRf*Err.e-F.KomegaF*omegaR;

xdot=[qrDot;omegaRDot];
end
