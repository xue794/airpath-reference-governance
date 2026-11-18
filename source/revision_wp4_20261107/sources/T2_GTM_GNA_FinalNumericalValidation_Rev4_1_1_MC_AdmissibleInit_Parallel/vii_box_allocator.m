function O=vii_box_allocator(B,tauReq,lo,hi,W,tol)
%VII_BOX_ALLOCATOR Exact active-set enumeration for 3-variable box LS.
%
% min_delta || W (B delta - tauReq) ||_2^2
% subject to lo <= delta <= hi.
%
% No Optimization Toolbox is required.  Every variable is either free,
% fixed at its lower bound, or fixed at its upper bound: 3^3=27 cases.

if nargin<6, tol=1e-11; end
if nargin<5 || isempty(W), W=eye(3); end

bestCost=inf;
bestNorm=inf;
bestDelta=zeros(3,1);
bestStatus=zeros(3,1);

states=[-1 0 1];

for i1=1:3
    for i2=1:3
        for i3=1:3
            st=[states(i1);states(i2);states(i3)];
            d=zeros(3,1);

            fixed=(st~=0);
            free=(st==0);

            d(st==-1)=lo(st==-1);
            d(st== 1)=hi(st== 1);

            rhs=tauReq-B(:,fixed)*d(fixed);

            if any(free)
                A=W*B(:,free);
                b=W*rhs;
                if rank(A)==sum(free)
                    df=A\b;
                else
                    df=pinv(A)*b;
                end
                d(free)=df;
            end

            if any(d<lo-tol) || any(d>hi+tol)
                continue;
            end

            d=min(max(d,lo),hi);
            r=B*d-tauReq;
            cost=norm(W*r)^2;
            dn=norm(d);

            if (cost<bestCost-1e-14) || ...
               (abs(cost-bestCost)<=1e-14 && dn<bestNorm)
                bestCost=cost;
                bestNorm=dn;
                bestDelta=d;
                bestStatus=st;
            end
        end
    end
end

if ~isfinite(bestCost)
    error('vii_box_allocator:noCandidate','No box-LS candidate found.');
end

O.delta=bestDelta;
O.residual=B*bestDelta-tauReq;
O.cost=bestCost;
O.status=bestStatus;
O.active=bestStatus~=0;
end
