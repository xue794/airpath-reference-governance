function C = t2_gna_coeff(alpha,beta,omega,delta,V,P)
%T2_GNA_COEFF  Grauer--Morelli 45-parameter GNA coefficient model.
% Inputs alpha,beta,delta in rad; omega=[wx wy wz]' in rad/s; V in m/s.
% Output C = [CD CY CL Cl Cm Cn]'.

if V <= 0
    error('t2_gna_coeff:V','Airspeed must be positive.');
end
p = omega(1); q = omega(2); r = omega(3);
da = delta(1); de = delta(2); dr = delta(3);
pt = P.b*p/(2*V);
qt = P.cbar*q/(2*V);
rt = P.b*r/(2*V);
t = P.theta;

CD = t(1)+t(2)*alpha+t(3)*alpha*qt+t(4)*alpha*de+t(5)*alpha^2 ...
   + t(6)*alpha^2*qt+t(7)*alpha^2*de+t(8)*alpha^3 ...
   + t(9)*alpha^3*qt+t(10)*alpha^4;
CY = t(11)*beta+t(12)*pt+t(13)*rt+t(14)*da+t(15)*dr;
CL = t(16)+t(17)*alpha+t(18)*qt+t(19)*de+t(20)*alpha*qt ...
   + t(21)*alpha^2+t(22)*alpha^3+t(23)*alpha^4;
Cl = t(24)*beta+t(25)*pt+t(26)*rt+t(27)*da+t(28)*dr;
Cm = t(29)+t(30)*alpha+t(31)*qt+t(32)*de+t(33)*alpha*qt ...
   + t(34)*alpha^2*qt+t(35)*alpha^2*de+t(36)*alpha^3*qt ...
   + t(37)*alpha^3*de+t(38)*alpha^4;
Cn = t(39)*beta+t(40)*pt+t(41)*rt+t(42)*da+t(43)*dr ...
   + t(44)*beta^2+t(45)*beta^3;
C = [CD;CY;CL;Cl;Cm;Cn];
end
