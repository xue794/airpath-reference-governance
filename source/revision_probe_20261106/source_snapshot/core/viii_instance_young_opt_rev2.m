function Y=viii_instance_young_opt_rev2(A,C,P8,P2)
%VIII_INSTANCE_YOUNG_OPT_REV2 Proof-weight optimization for Eq. (125)-(132).
%
% Controller gains are fixed.  Only Young inequality bookkeeping is varied.
% Self-dependence of each nonlinear interface is first absorbed into the
% receiving row's nominal dissipation.

L=A.L;
mu0=P2.mu0(:);
b=P2.b(:);

% Analytic dependency graph.
deps=cell(8,1);
deps{1}=[2 5];
deps{2}=[3 4];
deps{3}=[];
deps{4}=[6 7 8];
deps{5}=[1 2 3 4];
deps{6}=[1 2 3 4 5 7 8];
deps{7}=[1 2 3 4 5 6 8];
deps{8}=[1 2 3 4 5 6 7];

% Absorb sampled self-interface slopes for rows p,G,F,R,A.
self=zeros(8,1);
self(4)=L(4,4);
self(5)=L(5,5);
self(6)=L(6,6);
self(7)=L(7,7);
self(8)=L(8,8);

muNet=mu0-b.*self;

Y.mu0=mu0;
Y.selfL=self;
Y.muNet=muNet;
Y.validSelf=all(muNet>0);

if ~Y.validSelf
    Y.rho=inf;
    Y.D=diag(max(muNet,0));
    Y.G=nan(8);
    Y.retention=nan(8,1);
    Y.message='At least one sampled self-interface slope exhausts nominal dissipation.';
    return;
end

grid=P2.young.retentionGrid;
ext=P2.young.externalLossFraction;

% Initial retention fractions.
r=0.45*ones(8,1);
r(3)=1-ext; % wind row has no internal coupling.

% Coordinate descent over proof-only retention fractions.
for pass=1:P2.young.maxPasses
    changed=false;
    for i=1:8
        if isempty(deps{i})
            continue;
        end

        bestRi=r(i);
        bestRho=inf;

        for kk=1:numel(grid)
            rt=grid(kk);
            if rt>=1-ext
                continue;
            end
            rtry=r;
            rtry(i)=rt;
            [~,~,rho]=local_build(muNet,b,L,deps,rtry,ext);

            if isfinite(rho) && rho<bestRho
                bestRho=rho;
                bestRi=rt;
            end
        end

        if abs(bestRi-r(i))>1e-12
            changed=true;
            r(i)=bestRi;
        end
    end

    if ~changed
        break;
    end
end

[D,G,rho,epsCell]=local_build(muNet,b,L,deps,r,ext);

Y.D=D;
Y.G=G;
Y.rho=rho;
Y.retention=r;
Y.epsCell=epsCell;
Y.message='Young weights optimized with frozen controller parameters.';

if rho<1
    M=D-G;
    v=M\ones(8,1);
    w=M.'\ones(8,1);
    Y.M=M;
    Y.v=v;
    Y.w=w;
    Y.theta=w;
    Y.MmatrixPositive=all(v>0) && all(w>0);
else
    Y.M=nan(8);
    Y.v=nan(8,1);
    Y.w=nan(8,1);
    Y.theta=nan(8,1);
    Y.MmatrixPositive=false;
end
end

function [D,G,rho,epsCell]=local_build(muNet,b,L,deps,r,ext)
D=zeros(8);
G=zeros(8);
epsCell=cell(8,1);

for i=1:8
    if isempty(deps{i})
        D(i,i)=(1-ext)*muNet(i);
        epsCell{i}=[];
        continue;
    end

    ri=r(i);
    crossLoss=(1-ext-ri)*muNet(i);
    if crossLoss<=0
        rho=inf;
        D=nan(8); G=nan(8); epsCell={};
        return;
    end

    jj=deps{i};
    ell=L(i,jj);
    s=sum(ell);

    if s<=0
        D(i,i)=(1-ext)*muNet(i);
        epsCell{i}=zeros(size(ell));
        continue;
    end

    % sum eps_ij / 2 = crossLoss; row-optimal allocation for sum gamma.
    epsv=2*crossLoss*(ell/s);
    epsCell{i}=epsv;
    D(i,i)=ri*muNet(i);

    for k=1:numel(jj)
        if ell(k)>0
            G(i,jj(k))=(b(i)^2*ell(k)^2)/(2*epsv(k));
        end
    end
end

rho=max(abs(eig(D\G)));
end
