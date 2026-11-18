function xdot=viC_spin_rhs(t,x,F)
%VIC_SPIN_RHS Persistent pure-yaw command Cc=Rz(Omega t).

Cc=viB_Rz(F.spinRate*t);
qr=x(1:4);
omegaR=x(5:7);
Cr=quat_to_dcm(qr);
Q=Cc.'*Cr;
e=0.5*viC_vee(Q-Q.');

xdot=[quat_kinematics(qr,omegaR,F.quatNormGain);
      -F.kRf*e-F.KomegaF*omegaR];
end
