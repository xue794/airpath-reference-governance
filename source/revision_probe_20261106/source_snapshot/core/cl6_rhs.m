function xdot=cl6_rhs(t,x,mode,C)
%CL6_RHS Complete wind-aware T-2 GTM--GNA 6DOF closed-loop dynamics.

S=cl6_control_snapshot(t,x,mode,C);

% Reuse the already-verified Rev.1 quaternion plant directly.
U=cl6_plant_input(x(33:35),x(36),mode,C);
xpDot=t2_gna_6dof_quat(t,x(1:13),U,C.P);

xdot=[xpDot;
      S.zeDot;
      S.whatDot;
      S.xfDot;
      S.qrDot;
      S.omegaRDot;
      S.deltaDot;
      S.TDot];
end
