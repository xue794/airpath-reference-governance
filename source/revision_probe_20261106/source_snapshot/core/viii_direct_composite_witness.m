function W=viii_direct_composite_witness(Data,C,P8,J)
%VIII_DIRECT_COMPOSITE_WITNESS Direct full reduced-state sampled quadratic witness.
%
% This is a fallback numerical witness when no two-macro split is self-dissipative.
% It retains ALL signed internal couplings in the 34-state normalized Jacobian.

if nargin<4 || isempty(J)
    J=viii_collect_macro_jacobians(Data,C,P8);
end

G=viii_common_group_storage(J.A,P8);
aRaw=G.decayRaw;
a=P8.direct.selfDeflation*aRaw;

W.mode='direct';
W.valid=(a>P8.direct.minDecay);
W.storage=G;
W.a=a;
W.D=a;
W.G=0;
W.rho=0;
W.P=G.P;
W.J=J;
W.nWitnessPoints=numel(J.A);

if W.valid
    W.M=a;
    W.pi=1;
    W.Pi=1;
    W.Q=2*a;
    W.minEigQ=2*a;
else
    W.M=nan;
    W.pi=nan;
    W.Pi=nan;
    W.Q=nan;
    W.minEigQ=nan;
end
end
