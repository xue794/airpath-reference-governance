function O=viii_instance_interface_values_rev2(t,X,mode,C,P8,P2)
%VIII_INSTANCE_INTERFACE_VALUES_REV2 Interface signals used by Eq. (124a).
%
% The maps returned here correspond to the receiving rows:
%   p : d_p in Eq. (100)-(102)
%   G : col(-dot c_c, r_u) in Eq. (68a)-(68d)
%   F : omega_c in Eq. (95i)-(95k)
%   R : physical d_tau^B in Eq. (99e),(115),(119b)
%   A : dot u_A^star in Eq. (113c)-(113e)
%
% They are extracted from the SAME frozen nonlinear implementation, but the
% block graph itself remains the analytic graph from the theorem.

if nargin<6 || isempty(P2), P2=viii_instance_gain_params_rev2(C,P8); end

B=viii_blocks();
[fX,~]=viii_reduced_vector_field(t,X,mode,C,P8);
[x,~]=viii_lift_reduced(t,X,mode,C,P8);
[~,R0]=viii_reduced_from_full(t,x,mode,C,P8);
S=R0.S;

% p-row residual: xi_p_dot = -Kp xi_p + d_p.
xp=X(B.idx{4});
Kp=diag([C.A.kV C.A.kGamma C.A.kChi]);
dp=fX(B.idx{4})+Kp*xp;

% G-row affine perturbation: xG_dot = AG xG + uG.
xG=X(B.idx{5});
AG=[zeros(3),eye(3);-(C.R.Omega^2),-2*C.R.Zeta*C.R.Omega];
uG=fX(B.idx{5})-AG*xG;

% F-row command angular velocity.  From
% eRf_dot = E(Qf)(omega_r-Qf'*omega_c).
xF=X(B.idx{6});
Qf=viii_error_vector_to_dcm(xF(1:3));
E=0.5*(trace(Qf)*eye(3)-Qf.');
de=fX(B.idx{6}(1:3))-E*xF(4:6);
if rcond(E)<1e-10
    error('viii_instance_interface_values_rev2:E','E(Qf) became ill-conditioned.');
end
omegaC=-Qf*(E\de);

% R-row exact physical moment residual relative to the VI-D nominal target.
[~,tauA]=t2_gna_aero(S.Air.vgB,S.Air.omega,x(33:35), ...
    S.Air.CIB,S.Air.wI,S.Air.h,C.P);
dTau=tauA-(-S.Mom.D*S.Air.omega+S.tauC);

% A-row command rate.  xA_dot = actuator_dot - uAstar_dot.
actDot=[S.deltaDot;S.TDot];
uADot=actDot-fX(B.idx{8});

O.p=dp;
O.G=uG;
O.F=omegaC;
O.R=dTau;
O.A=uADot;
O.fX=fX;
O.S=S;
end
