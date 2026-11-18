function P3=viii_instance_gain_params_rev3(C,P8)
%VIII_INSTANCE_GAIN_PARAMS_REV3 Structured self-dissipation instance settings.
%
% Paper logic:
%   the repaired theorem is primary;
%   this file only instantiates the theorem for the frozen T-2 GTM--GNA case.
%
% No controller, plant, actuator, or theorem parameter is retuned.

if nargin<1 || isempty(C), C=cl6_params(); end
if nargin<2 || isempty(P8), P8=viii_params(); end

P3.revision='Analytic Gain Ledger Numerical Instance Verification Rev.3';

P3.relStep=2e-4;
P3.maxHalves=10;
P3.interfaceInflation=1.10;

P3.sample.hold=[10 40];
P3.sample.varying=[15 35 55 75 95 115];
P3.sample.noisy=[20 40 60 80 100 115];

% Reserve a small fraction of each sampled self-dissipation margin for
% primitive/reference forcing before the internal Young split.
P3.young.externalLossFraction=0.01;
P3.young.retentionGrid=0.05:0.05:0.90;
P3.young.maxPasses=10;

% Governor storage.
AG=[zeros(3),eye(3);-(C.R.Omega^2),-2*C.R.Zeta*C.R.Omega];
PG=rg_solve_lyap(AG,eye(6));
P3.AG=AG;
P3.PG=PG;
P3.QG=eye(6);

% Exact direct z/s/w rows.
P3.alphaExact=zeros(8,1);
P3.alphaExact(1)=min(diag(C.G.Lambda));
P3.alphaExact(2)=min(diag(C.G.Ks));
P3.alphaExact(3)=min(eig(C.E.Lw));

P3.cExact=zeros(8,8);
P3.cExact(1,2)=1;
P3.cExact(1,5)=norm(C.G.Lambda\eye(3),2)*C.R.Vmax;

Vstar=P8.Vstar;
Lvp=max([1,C.R.Vmax/Vstar]);
P3.cExact(2,3)=1;
P3.cExact(2,4)=Lvp;

P3.labels={'z','s','w','p','G','F','R','A'};
end
