function G=viB_airdata_from_geometry(C,t)
%VIB_AIRDATA_FROM_GEOMETRY Extract alpha,beta from inertial air-path direction.
%
% b=C'*t is the unit air-relative velocity direction in actual body axes:
%   b=[cos(alpha)cos(beta); sin(beta); sin(alpha)cos(beta)].

b=C.'*t;
b=b/norm(b);
s=max(-1,min(1,b(2)));
beta=asin(s);
alpha=atan2(b(3),b(1));

G.b=b;
G.alpha=alpha;
G.beta=beta;
G.kappa=cos(alpha)*cos(beta);
end
