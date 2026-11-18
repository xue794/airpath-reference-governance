function U=viii_macro_forcing_uub(Tv,Tn,W,P8)
%VIII_MACRO_FORCING_UUB Empirical forcing/remainder envelope and UUB comparison.
%
% The formal theorem requires a continuous compact-domain forcing/remainder
% envelope.  Numerically, Rev.1.1 estimates the required envelope directly from
% the accepted varying/noisy nonlinear trajectories:
%
%   h_req = max(0, rdot + D r - G r).
%
% The maximum is inflated by the frozen witness safety factor.

a=W.a(:);
G=W.c;

Tset={Tv,Tn};
hmax=zeros(2,1);

for it=1:2
    T=Tset{it};
    for k=1:size(T.r,1)
        rv=T.r(k,:).';
        rd=T.rdot(k,:).';
        hreq=rd+diag(a)*rv-G*rv;
        hmax=max(hmax,max(hreq,0));
    end
end

hbar=P8.macro.forcingInflation*hmax;
M=diag(a)-G;

rInf=M\hbar;

% Finite-horizon positive comparison from each actual initial root state.
sets={Tv,Tn};
maxRatio=0;
maxViolation=0;
Bnd=cell(2,1);

for it=1:2
    T=sets{it};
    n=size(T.r,1);
    rb=zeros(n,2);
    r0=T.r(1,:).';

    for k=1:n
        tt=Tset_time(T,k);
        E=expm(-M*tt);
        rb(k,:)=(E*r0 + M\(eye(2)-E)*hbar).';
    end

    Bnd{it}=rb;
    ratio=T.r./max(rb,1e-12);
    maxRatio=max(maxRatio,max(ratio(:)));
    tmp=T.r-rb;
    maxViolation=max(maxViolation,max(tmp(:)));
end

U.hRequired=hmax;
U.hbar=hbar;
U.rInf=rInf;
U.bounds=Bnd;
U.maxRatio=maxRatio;
U.maxViolation=maxViolation;
end

function tt=Tset_time(T,k)
% T itself stores no time vector; recover time from the uniform sample index
% using a private field if present, otherwise assume caller added it.
if isfield(T,'t')
    tt=T.t(k)-T.t(1);
else
    error('viii_macro_forcing_uub:Time','Trajectory storage must contain time.');
end
end
