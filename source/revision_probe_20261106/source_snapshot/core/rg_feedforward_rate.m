function D = rg_feedforward_rate(x,u,R)
%RG_FEEDFORWARD_RATE Analytic first derivatives; affine in u.
V=x(1); gamma=x(2); nuV=x(4); nuG=x(5); nuC=x(6);
F=rg_feedforward(x,R);
if ~F.rootInfo.bracket
    error('rg_feedforward_rate:root','Feedforward AoA root not bracketed.');
end

dag = nuV*nuG + V*u(2) - R.g*sin(gamma)*nuG;
dac = nuV*cos(gamma)*nuC - V*sin(gamma)*nuG*nuC + V*cos(gamma)*u(3);
Xdot = R.P.m*(u(1)+R.g*cos(gamma)*nuG);

rr=hypot(F.aGamma,F.aChi);
Ndot=R.P.m*(F.aGamma*dag+F.aChi*dac)/rr;

M=rg_N_map(F.alpha,V,F.X,R);
alphaDot=(Ndot-M.V*nuV-M.X*Xdot)/M.alpha;

A=rg_nominal_aero(V,F.alpha,R);
Ddot=A.DV*nuV+A.Da*alphaDot;
Tdot=(Ddot+Xdot)/cos(F.alpha) + ...
     (A.D+F.X)/cos(F.alpha)*tan(F.alpha)*alphaDot;

D.dag=dag; D.dac=dac; D.Xdot=Xdot; D.Ndot=Ndot;
D.alphaDot=alphaDot; D.Tdot=Tdot;
end
