function S=viii_structured_storage_coeff_rev3(X,O,C,P8,P3)
%VIII_STRUCTURED_STORAGE_COEFF_REV3 Exact receiving-storage coefficients.
%
% For i = p,G,F,R,A this returns matrices satisfying the local storage form
%
%   Vdot_i <= -x_i' R_i x_i + x_i' N_i u_i + primitive/reference terms.
%
% If u_i = H_ii x_i + sum_j H_ij x_j + d_i, then
%
%   Q_i,net = R_i - sym(N_i H_ii)
%
% is the STRUCTURED self-dissipation matrix.  This is the quantity that must
% be checked before Eq. (124a); it is not replaced by
%   lambda_min(R_i) - ||N_i|| ||H_ii||.

B=viii_blocks();

S.R=cell(8,1);
S.N=cell(8,1);

% p row: Vp=0.5||xi_p||^2,
% xi_p_dot=-Kp xi_p+d_p.
Kp=diag([C.A.kV C.A.kGamma C.A.kChi]);
S.R{4}=Kp;
S.N{4}=eye(3);

% G row: VG=xG'PG xG and AG'PG+PG AG=-I.
S.R{5}=P3.QG;
S.N{5}=2*P3.PG;

% F row: exact SO(3) storage (95d).
xF=X(B.idx{6});
eF=xF(1:3);
wR=xF(4:6);
Qf=viii_error_vector_to_dcm(eF);
Ef=0.5*(trace(Qf)*eye(3)-Qf.');

KF=C.F.KomegaF;
cF=C.F.cF;
kF=C.F.kRf;

RF=[cF*kF*eye(3),          (cF/2)*KF; ...
    (cF/2)*KF, KF-cF*0.5*(Ef+Ef.')];

% omega_c input coefficient:
% Vdot_F,input = -[kF*Qf*eF + cF*Qf*Ef'*wR]' omega_c
%               = xF' NF omega_c.
BF=[kF*Qf, cF*Qf*Ef.'];
NF=-BF.';

S.R{6}=0.5*(RF+RF.');
S.N{6}=NF;

% R row: exact rotational storage (99a).
xR=X(B.idx{7});
eR=xR(1:3);
QR=viii_error_vector_to_dcm(eR);
ER=0.5*(trace(QR)*eye(3)-QR.');

Kd=O.S.Ctl.Kd;
I=C.P.I;
cR=C.Dp.cR;
kR=C.Dp.kR;

RR=[cR*kR*eye(3),           (cR/2)*Kd; ...
    (cR/2)*Kd, Kd-cR*0.5*(ER.'*I+I*ER)];

% Vdot_R,input = (eOmega+cR eR)' d_tau = xR' NR d_tau.
CR=[cR*eye(3),eye(3)];
NR=CR.';

S.R{7}=0.5*(RR+RR.');
S.N{7}=NR;

% A row: VA=0.5||xA||^2 and the rate-saturated projected actuator
% sector gives xA' actuator_restoring <= -xA' KappaA xA.
kap=[P8.kappaDelta(:);P8.kappaT];
S.R{8}=diag(kap);

% xA_dot = actuator_restoring - dot(u_A^star),
% so Vdot_A,input = -xA' dot(u_A^star).
S.N{8}=-eye(4);
end
