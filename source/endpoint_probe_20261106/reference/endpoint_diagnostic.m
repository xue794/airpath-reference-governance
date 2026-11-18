function d=endpoint_diagnostic(x,R)
F=rg_feedforward(x,R);V=x(1);ga=x(2);nu=x(4:6);
raw=[F.hNminus+R.DeltaN;F.hNplus+R.DeltaN];
psi=[nu(1)+R.k1*(V-R.Vmin);-nu(1)+R.k1*(R.Vmax-V); ...
nu(2)+R.k1*(ga+R.gammaMax);-nu(2)+R.k1*(R.gammaMax-ga)];
con=[V-R.Vmin;R.Vmax-V;ga+R.gammaMax;R.gammaMax-ga;R.nuMax-abs(nu); ...
F.aGamma-R.aGammaMin;R.muMax-abs(F.mu);psi];
% Endpoint values of thrust are well-defined even at a root event.
Al=rg_nominal_aero(V,R.alphaLo,R);Ah=rg_nominal_aero(V,R.alphaHi,R);
Tl=(Al.D+F.X)/cos(R.alphaLo);Th=(Ah.D+F.X)/cos(R.alphaHi);
d=struct('raw',raw,'reserve',raw-10,'T',F.T,'root',F.rootInfo.bracket,'alpha',F.alpha, ...
'conventionalMargins',con,'conventionalOK',all(con>=-1e-8),'psi',psi, ...
'endpointThrust',[Tl;Th],'bank',F.mu);
end
