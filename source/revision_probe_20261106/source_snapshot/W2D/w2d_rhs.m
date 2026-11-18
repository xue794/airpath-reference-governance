function xdot=w2d_rhs(t,x,mode,C,method,O)
%W2D_RHS Full 6DOF dynamics; all targeted cases use the source still-wind plant mode.
S=w2d_control_snapshot(t,x,mode,C,method,O);U=cl6_plant_input(x(33:35),x(36),'hold',C);xpDot=t2_gna_6dof_quat(t,x(1:13),U,C.P);
xdot=[xpDot;S.zeDot;S.whatDot;S.xfDot;S.qrDot;S.omegaRDot;S.deltaDot;S.TDot];
end
