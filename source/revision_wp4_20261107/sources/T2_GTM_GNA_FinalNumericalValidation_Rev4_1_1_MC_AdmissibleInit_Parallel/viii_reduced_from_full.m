function [X,O]=viii_reduced_from_full(t,x,mode,C,P8)
%VIII_REDUCED_FROM_FULL Compute the 34-state reduced tracking coordinate.
%
% X = col(z_e,s,wtilde,xi_p,x_G,x_F,x_R,x_A)

if nargin<5 || isempty(P8), P8=viii_params(); end

S=cl6_control_snapshot(t,x,mode,C);
xf=x(20:25);
omegaR=x(30:32);

z=x(14:16);
s=S.Gd.s;
wt=S.tildew;

eV=S.Air.V-xf(1);
egamma=S.Air.gamma-xf(2);
echi=rg_angle_error(S.Air.chi,xf(3));
xp=[eV;P8.Vstar*egamma;P8.Vstar*echi];

eta=[xf(1)-S.cc(1);
     xf(2)-S.cc(2);
     rg_angle_error(xf(3),S.cc(3))];
xG=[eta;xf(4:6)];

xF=[S.eRf;omegaR];

ER=viD_error(S.Cr,omegaR,S.Air.CIB,S.Air.omega);
xR=[ER.eR;ER.eOmega];

xA=[x(33:35)-S.deltaCmd;
    x(36)-S.Tcmd];

X=[z;s;wt;xp;xG;xF;xR;xA];

O.S=S;
O.ER=ER;
O.blocks=viii_blocks();
end
