function M=vii_moment_components(V,alpha,beta,R)
%VII_MOMENT_COMPONENTS Exact nominal T-2 GNA moment decomposition.
%
% tau_a = tau0 + (-Domega+Somega)*omega + Bdelta*delta.

th=R.P.theta;
qS=0.5*R.rho*R.P.S*V^2;

Cl0=th(24)*beta;
Cm0=th(29)+th(30)*alpha+th(38)*alpha^4;
Cn0=th(39)*beta+th(44)*beta^2+th(45)*beta^3;

tau0=qS*[R.P.b*Cl0;
         R.P.cbar*Cm0;
         R.P.b*Cn0];

Rate=viD_rate_moment_decomp(V,alpha,R);
Bdelta=vii_surface_effectiveness(V,alpha,R);

M.tau0=tau0;
M.K=Rate.K;
M.D=Rate.D;
M.S=Rate.S;
M.B=Bdelta;
end
