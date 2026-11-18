function fDelta=vii_surface_force_residual(V,alpha,delta,R)
%VII_SURFACE_FORCE_RESIDUAL Translational aerodynamic force caused by surfaces.
%
% This force is deliberately NOT cancelled by the rotational allocator.
% It belongs to the translational model-residual channel d_F,m in the paper.

th=R.P.theta;
qS=0.5*R.rho*R.P.S*V^2;

da=delta(1); de=delta(2); dr=delta(3);

dCD=(th(4)*alpha+th(7)*alpha^2)*de;
dCY=th(14)*da+th(15)*dr;
dCL=th(19)*de;

fDelta=qS*[ -dCD*cos(alpha)+dCL*sin(alpha);
             dCY;
            -dCD*sin(alpha)-dCL*cos(alpha) ];
end
