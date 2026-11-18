function [fB,tauB,aux] = t2_gna_aero(vgB,omega,delta,CIB,wI,h,P)
%T2_GNA_AERO  Wind-aware GNA aerodynamic forces/moments in study body axes.
% Body axes are x forward, y right, z down. Positive CL is implemented as
% upward aerodynamic lift, hence negative body-z force near alpha=0.

vaB = vgB - CIB.'*wI;
V = norm(vaB);
if V < 1e-6
    error('t2_gna_aero:V','Airspeed too small for GNA model.');
end
alpha = atan2(vaB(3),vaB(1));
beta = asin(max(-1,min(1,vaB(2)/V)));
rho = t2_isa_density(h);
qbar = 0.5*rho*V^2;
C = t2_gna_coeff(alpha,beta,omega,delta,V,P);
CD=C(1); CY=C(2); CL=C(3); Cl=C(4); Cm=C(5); Cn=C(6);

% Force mapping adopted for the paper's x-forward/y-right/z-down body axes.
% It preserves the standard physical signs: drag opposes +x at small alpha,
% positive lift acts in -z, and CY acts in +y.
fB = qbar*P.S*[ -CD*cos(alpha)+CL*sin(alpha); ...
                 CY; ...
                -CD*sin(alpha)-CL*cos(alpha) ];
tauB = qbar*P.S*[P.b*Cl; P.cbar*Cm; P.b*Cn];

aux = struct('vaB',vaB,'V',V,'alpha',alpha,'beta',beta,'rho',rho, ...
             'qbar',qbar,'C',C);
end
