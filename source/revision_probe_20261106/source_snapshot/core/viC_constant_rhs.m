function xdot=viC_constant_rhs(~,x,Cc,F)
%VIC_CONSTANT_RHS Exact constant-command filter simulation.
qr=x(1:4);
omegaR=x(5:7);
Cr=quat_to_dcm(qr);
Q=Cc.'*Cr;
e=0.5*viC_vee(Q-Q.');

xdot=[quat_kinematics(qr,omegaR,F.quatNormGain);
      -F.kRf*e-F.KomegaF*omegaR];
end
