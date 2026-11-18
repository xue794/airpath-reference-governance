function xdot=w2c_rhs(t,x,mode,C,method,O)
%W2C_RHS Full 6DOF dynamics with the W2-C governor-selection interface.
S=w2c_control_snapshot(t,x,mode,C,method,O);
U=cl6_plant_input(x(33:35),x(36),mode,C);
xpDot=t2_gna_6dof_quat(t,x(1:13),U,C.P);
xdot=[xpDot;S.zeDot;S.whatDot;S.xfDot;S.qrDot;S.omegaRDot;S.deltaDot;S.TDot];
end
