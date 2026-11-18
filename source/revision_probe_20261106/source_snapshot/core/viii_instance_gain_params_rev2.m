function P2=viii_instance_gain_params_rev2(C,P8)
%VIII_INSTANCE_GAIN_PARAMS_REV2 Numerical instance-verification settings.
%
% This file does NOT modify the theorem and does NOT tune the controller.
% It instantiates the analytic gain ledger on feasible points of the frozen
% T-2 GTM--GNA study case.
%
% The numerical audit estimates local blockwise Lipschitz envelopes of the
% explicit analytic interfaces used in Eq. (124a)-(132).  These estimates
% are a numerical instance check, not a replacement for the theorem.

if nargin<1 || isempty(C), C=cl6_params(); end
if nargin<2 || isempty(P8), P8=viii_params(); end

P2.revision='Analytic Gain Ledger Numerical Instance Verification Rev.2';

% Central block perturbation.  Actual coordinate step is relStep times the
% physical scale already used by the Section-VIII reduced-coordinate tools.
P2.relStep=2e-4;
P2.maxHalves=10;

% Multiplicative margin applied AFTER taking the maximum over all feasible
% centers.  This is a numerical-validation margin, not a theorem constant.
P2.interfaceInflation=1.10;

% Only a moderate set of feasible centers is needed for the instance audit.
% Expensive full trajectories are cached / reused.
P2.sample.hold=[10 40];
P2.sample.varying=[15 35 55 75 95 115];
P2.sample.noisy=[20 40 60 80 100 115];

% Young-weight search.  This changes proof bookkeeping only.
P2.young.externalLossFraction=0.01;
P2.young.retentionGrid=0.05:0.05:0.90;
P2.young.maxPasses=8;

% Direct analytic constants.
P2.mu0=zeros(8,1);
P2.b=zeros(8,1);

P2.mu0(1)=min(diag(C.G.Lambda));
P2.b(1)=1;

P2.mu0(2)=min(diag(C.G.Ks));
P2.b(2)=1;

P2.mu0(3)=min(eig(C.E.Lw));
P2.b(3)=1;

P2.mu0(4)=min([C.A.kV C.A.kGamma C.A.kChi]);
P2.b(4)=1;

AG=[zeros(3),eye(3);-(C.R.Omega^2),-2*C.R.Zeta*C.R.Omega];
PG=rg_solve_lyap(AG,eye(6));
P2.PG=PG;
P2.mu0(5)=1;
P2.b(5)=2*norm(PG,2);

kbarF=norm(C.F.KomegaF,2);
kminF=min(eig(C.F.KomegaF));
MF=[C.F.cF*C.F.kRf, -C.F.cF*kbarF/2; ...
   -C.F.cF*kbarF/2, kminF-C.F.cF];
P2.MF=MF;
P2.mu0(6)=min(eig(MF));
P2.b(6)=hypot(C.F.kRf,C.F.cF);

DC=viD_damping_certificate(C.R,C.Dp);
kminR=min(eig(C.Dp.KomegaC))+DC.minEigD;
kmaxR=norm(C.Dp.KomegaC,2)+DC.maxNormD;
Ibar=max(eig(C.P.I));
MR=[C.Dp.cR*C.Dp.kR, -C.Dp.cR*kmaxR/2; ...
   -C.Dp.cR*kmaxR/2, kminR-C.Dp.cR*Ibar];
P2.MR=MR;
P2.mu0(7)=min(eig(MR));
P2.b(7)=hypot(1,C.Dp.cR);

P2.mu0(8)=P8.kappaA;
P2.b(8)=1;

% Exact z/s analytic interfaces.
P2.Lexact=zeros(8,8);
P2.Lexact(1,2)=1;
P2.Lexact(1,5)=norm(C.G.Lambda\eye(3),2)*C.R.Vmax;

Vstar=P8.Vstar;
Lvp=max([1,C.R.Vmax/Vstar,C.R.Vmax/(Vstar)]); % cos(gamma)<=1
P2.Lexact(2,3)=1;
P2.Lexact(2,4)=Lvp;

P2.labels={'z','s','w','p','G','F','R','A'};
end
