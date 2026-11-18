function edot = viA_error_rhs_theory(y,xf,wdot,R,A)
%VIA_ERROR_RHS_THEORY Reduced error dynamics predicted by Section VI-A.
C=viA_feedback_command(y,xf,R,A);
B=viA_airpath_basis(y(2),y(3));
V=y(1); gamma=y(2);

edot=[-A.kV*C.e(1)-B.t.'*wdot;
      -A.kGamma*C.e(2)-B.nGamma.'*wdot/V;
      -A.kChi*C.e(3)-B.nChi.'*wdot/(V*cos(gamma))];
end
