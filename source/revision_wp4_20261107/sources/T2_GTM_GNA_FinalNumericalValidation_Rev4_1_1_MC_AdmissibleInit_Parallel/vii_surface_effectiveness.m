function B=vii_surface_effectiveness(V,alpha,R)
%VII_SURFACE_EFFECTIVENESS Exact T-2 GNA moment effectiveness wrt [da de dr].
%
% tau_delta = B_delta(V,alpha)*delta.

th=R.P.theta;
qS=0.5*R.rho*R.P.S*V^2;

B=zeros(3,3);

B(1,1)=qS*R.P.b*th(27);
B(1,3)=qS*R.P.b*th(28);

cmde=th(32)+th(35)*alpha^2+th(37)*alpha^3;
B(2,2)=qS*R.P.cbar*cmde;

B(3,1)=qS*R.P.b*th(42);
B(3,3)=qS*R.P.b*th(43);
end
