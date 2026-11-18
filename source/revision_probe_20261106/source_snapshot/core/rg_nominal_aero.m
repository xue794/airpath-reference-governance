function A = rg_nominal_aero(V,alpha,R)
%RG_NOMINAL_AERO D0,L0 and alpha derivatives for zero-rate, zero-surface GNA map.
th=R.P.theta;
CD = th(1)+th(2)*alpha+th(5)*alpha^2+th(8)*alpha^3+th(10)*alpha^4;
CL = th(16)+th(17)*alpha+th(21)*alpha^2+th(22)*alpha^3+th(23)*alpha^4;
CDa = th(2)+2*th(5)*alpha+3*th(8)*alpha^2+4*th(10)*alpha^3;
CLa = th(17)+2*th(21)*alpha+3*th(22)*alpha^2+4*th(23)*alpha^3;

qS = 0.5*R.rho*R.P.S*V^2;
A.CD=CD; A.CL=CL;
A.D=qS*CD; A.L=qS*CL;
A.Da=qS*CDa; A.La=qS*CLa;
A.DV=2*A.D/V; A.LV=2*A.L/V;
end
