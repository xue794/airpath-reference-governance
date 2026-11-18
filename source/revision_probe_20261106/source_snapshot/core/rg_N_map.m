function M = rg_N_map(alpha,V,X,R)
%RG_N_MAP Combined nominal lift/thrust normal-force map and derivatives.
A=rg_nominal_aero(V,alpha,R);
ta=tan(alpha);
M.value = A.L + (A.D+X)*ta;
M.alpha = A.La + A.Da*ta + (A.D+X)/(cos(alpha)^2);
M.V = A.LV + A.DV*ta;
M.X = ta;
M.aero=A;
end
