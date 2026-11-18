function xdot=viC_filter_interp_rhs(t,x,Cmd,F)
%VIC_FILTER_INTERP_RHS Command filter driven by interpolated VI-A attitude command.

Cc=viC_interp_command(t,Cmd);
qr=x(1:4);
omegaR=x(5:7);
Cr=quat_to_dcm(qr);

Q=Cc.'*Cr;
e=0.5*viC_vee(Q-Q.');

qrDot=quat_kinematics(qr,omegaR,F.quatNormGain);
omegaRDot=-F.kRf*e-F.KomegaF*omegaR;

xdot=[qrDot;omegaRDot];
end
