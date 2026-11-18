function U=viii_certificate_forcing_uub(Tv,Tn,W,P8)
%VIII_CERTIFICATE_FORCING_UUB Empirical forcing/remainder envelope.

if strcmp(W.mode,'macro')
    U=viii_macro_forcing_uub(Tv,Tn,W,P8);
    return;
end

a=W.a;
hmax=0;
Tset={Tv,Tn};

for it=1:2
    T=Tset{it};
    hreq=T.rdot+a*T.r;
    hmax=max(hmax,max(max(hreq,0)));
end

hbar=P8.direct.forcingInflation*hmax;
rInf=hbar/a;

Bnd=cell(2,1);
maxRatio=0;
maxViolation=0;

for it=1:2
    T=Tset{it};
    n=numel(T.t);
    rb=zeros(n,1);
    r0=T.r(1);

    for k=1:n
        tt=T.t(k)-T.t(1);
        rb(k)=exp(-a*tt)*r0+(1-exp(-a*tt))*rInf;
    end

    Bnd{it}=rb;
    ratio=T.r./max(rb,1e-12);
    maxRatio=max(maxRatio,max(ratio));
    maxViolation=max(maxViolation,max(T.r-rb));
end

U.hRequired=hmax;
U.hbar=hbar;
U.rInf=rInf;
U.bounds=Bnd;
U.maxRatio=maxRatio;
U.maxViolation=maxViolation;
end
