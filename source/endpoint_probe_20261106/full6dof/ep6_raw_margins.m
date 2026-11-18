function h=ep6_raw_margins(t,x,C)
% Four root-free quantities for TERMINAL events, original 0..9deg interval.
R=C.R;xf=x(20:25);V=xf(1);gam=xf(2);
X=R.P.m*(xf(4)+R.g*sin(gam));ag=V*xf(5)+R.g*cos(gam);ac=V*cos(gam)*xf(6);
N=R.P.m*hypot(ag,ac);lo=rg_N_map(R.alphaLo,V,X,R);hi=rg_N_map(R.alphaHi,V,X,R);
Air=cl6_airpath_state(t,x(1:13),'hold',C);A=C.A;
Xc=R.P.m*(xf(4)-A.kV*(Air.V-xf(1))+R.g*sin(Air.gamma));
ag=Air.V*(xf(5)-A.kGamma*(Air.gamma-xf(2)))+R.g*cos(Air.gamma);
ac=Air.V*cos(Air.gamma)*(xf(6)-A.kChi*rg_angle_error(Air.chi,xf(3)));
Nc=R.P.m*hypot(ag,ac);loC=rg_N_map(R.alphaLo,Air.V,Xc,R);hiC=rg_N_map(R.alphaHi,Air.V,Xc,R);
h=[N-lo.value;hi.value-N;Nc-loC.value;hiC.value-Nc];
end
