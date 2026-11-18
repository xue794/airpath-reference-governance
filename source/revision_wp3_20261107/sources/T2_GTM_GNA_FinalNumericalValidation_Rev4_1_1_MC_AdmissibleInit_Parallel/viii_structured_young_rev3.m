function Y=viii_structured_young_rev3(A,P3)
%VIII_STRUCTURED_YOUNG_REV3 Eq. (124a)->(125)->(132) instance bookkeeping.
%
% Inputs:
%   A.alpha(i) = sampled structured self-dissipation margin alpha_i
%   A.c(i,j)   = sampled bilinear cross coefficient c_ij in
%                c_ij ||x_i|| ||x_j||
%
% Only Young weights are optimized; controller gains are frozen.

alpha=A.alpha(:);
C=A.c;

deps=cell(8,1);
deps{1}=[2 5];
deps{2}=[3 4];
deps{3}=[];
deps{4}=[6 7 8];
deps{5}=[1 2 3 4];
deps{6}=[1 2 3 4 5 7 8];
deps{7}=[1 2 3 4 5 6 8];
deps{8}=[1 2 3 4 5 6 7];

Y.validSelf=all(alpha>0);
Y.alpha=alpha;

if ~Y.validSelf
    Y.rho=inf;
    Y.D=diag(max(alpha,0));
    Y.G=nan(8);
    Y.retention=nan(8,1);
    Y.message='At least one STRUCTURED sampled self-dissipation margin is nonpositive.';
    return;
end

ext=P3.young.externalLossFraction;
grid=P3.young.retentionGrid;

r=0.45*ones(8,1);
r(3)=1-ext;

for pass=1:P3.young.maxPasses
    changed=false;

    for i=1:8
        if isempty(deps{i}), continue; end

        bestR=r(i);
        bestRho=inf;

        for kk=1:numel(grid)
            rt=grid(kk);
            if rt>=1-ext, continue; end

            rr=r;
            rr(i)=rt;
            [~,~,rho]=local_build(alpha,C,deps,rr,ext);

            if isfinite(rho) && rho<bestRho
                bestRho=rho;
                bestR=rt;
            end
        end

        if abs(bestR-r(i))>1e-12
            r(i)=bestR;
            changed=true;
        end
    end

    if ~changed, break; end
end

[D,G,rho,epsCell]=local_build(alpha,C,deps,r,ext);

Y.D=D;
Y.G=G;
Y.rho=rho;
Y.retention=r;
Y.epsCell=epsCell;
Y.message='Structured self absorption + cross-only Young ledger.';

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

function [D,G,rho,epsCell]=local_build(alpha,C,deps,r,ext)
D=zeros(8);
G=zeros(8);
epsCell=cell(8,1);

for i=1:8
    if isempty(deps{i})
        D(i,i)=(1-ext)*alpha(i);
        epsCell{i}=[];
        continue;
    end

    crossLoss=(1-ext-r(i))*alpha(i);
    if crossLoss<=0
        D=nan(8); G=nan(8); rho=inf; epsCell={};
        return;
    end

    jj=deps{i};
    cc=C(i,jj);
    s=sum(cc);

    D(i,i)=r(i)*alpha(i);

    if s<=0
        epsCell{i}=zeros(size(cc));
        continue;
    end

    % Under sum eps_ij/2 = crossLoss, eps_ij proportional to c_ij
    % minimizes the row sum of Young-generated squared-state gains.
    epsv=2*crossLoss*(cc/s);
    epsCell{i}=epsv;

    for k=1:numel(jj)
        if cc(k)>0
            G(i,jj(k))=cc(k)^2/(2*epsv(k));
        end
    end
end

rho=max(abs(eig(D\G)));
end
