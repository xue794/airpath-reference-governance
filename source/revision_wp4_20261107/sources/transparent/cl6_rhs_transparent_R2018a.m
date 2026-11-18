function xdot = cl6_rhs_transparent_R2018a(t,x,mode,C)
%CL6_RHS_TRANSPARENT_R2018A
% Same full closed-loop dynamics as cl6_rhs, but with transparent governor.

S = cl6_control_snapshot_transparent_R2018a(t,x,mode,C);

U = cl6_plant_input(x(33:35),x(36),mode,C);
xpDot = t2_gna_6dof_quat(t,x(1:13),U,C.P);

xdot = [xpDot;
        S.zeDot;
        S.whatDot;
        S.xfDot;
        S.qrDot;
        S.omegaRDot;
        S.deltaDot;
        S.TDot];
end
