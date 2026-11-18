function F = rg_feedforward(x,R)
%RG_FEEDFORWARD Reference force geometry and nominal AoA/thrust inversion.
% x=[V;gamma;chi;nuV;nuGamma;nuChi]
V=x(1); gamma=x(2); nuV=x(4); nuG=x(5); nuC=x(6);

F.X = R.P.m*(nuV + R.g*sin(gamma));
F.aGamma = V*nuG + R.g*cos(gamma);
F.aChi = V*cos(gamma)*nuC;
F.N = R.P.m*hypot(F.aGamma,F.aChi);

[F.alpha,F.rootInfo]=rg_alpha_root(V,F.X,F.N,R);
if F.rootInfo.bracket
    A=rg_nominal_aero(V,F.alpha,R);
    F.Fparallel=A.D+F.X;
    F.T=F.Fparallel/cos(F.alpha);
else
    F.Fparallel=NaN; F.T=NaN;
end
F.mu=atan2(F.aChi,F.aGamma);

Mlo=rg_N_map(R.alphaLo,V,F.X,R);
Mhi=rg_N_map(R.alphaHi,V,F.X,R);
F.hNminus=F.N-Mlo.value-R.DeltaN;
F.hNplus=Mhi.value-F.N-R.DeltaN;
end
