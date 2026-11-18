function W=viii_small_gain_witness(FR,C,P8)
%VIII_SMALL_GAIN_WITNESS Construct the root-energy M-matrix certificate.
%
% Numerical local witness in normalized reduced coordinates:
%
%   ydot = A(y,t)y + b(y,t) + nonlinear remainder.
%
% At each FEASIBLE center y_k, the local affine forcing is defined exactly as
%
%   b_k = f(y_k,t_k) - A_k y_k.
%
% The sampled subsystem inequalities are
%
%   Vdot_i <= -a0_i||y_i||^2
%             + sum_j c_ij||y_i||||y_j||
%             + b_i||y_i||.
%
% With r_i=sqrt(V_i),
%
%   D^+ r <= -(D-G)r + h.
%
% The formal manuscript theorem remains condition-based on continuous
% compact-domain bounds; this function is a deterministic trajectory-tube
% numerical witness.

SM=viii_storage_matrices(C,P8);
B=SM.blocks;
mask=viii_dependency_mask();
Pts=viii_collect_witness_points(FR,C,P8);

a0=inf(8,1);
c=zeros(8,8);
maxForbidden=0;
hCenterMax=zeros(8,1);
roundTripMax=0;

for kp=1:numel(Pts)
    pt=Pts(kp);

    % Lifting round-trip audit only at FEASIBLE centers.
    [xLift,~]=viii_lift_reduced(pt.t,pt.X,pt.mode,C,P8);
    Xback=viii_reduced_from_full(pt.t,xLift,pt.mode,C,P8);
    roundTripMax=max(roundTripMax, ...
        norm((Xback-pt.X)./P8.scaleVector,inf));

    J=viii_local_jacobian(pt.t,pt.X,pt.mode,C,P8);
    A=J.A;

    % Feasible-center affine/reference forcing in normalized coordinates.
    bAff=J.f-A*J.y;

    for i=1:8
        ii=B.idx{i};
        Aii=A(ii,ii);
        Pi=SM.P{i};

        Qi=-0.5*(Pi*Aii+Aii.'*Pi);
        a0(i)=min(a0(i),min(eig(0.5*(Qi+Qi.'))));

        for j=1:8
            if j==i, continue; end
            jj=B.idx{j};
            cij=norm(Pi*A(ii,jj),2);

            if mask(i,j)
                c(i,j)=max(c(i,j),cij);
            else
                maxForbidden=max(maxForbidden,cij);
            end
        end

        hi=norm(Pi*bAff(ii))/(2*sqrt(SM.m(i)));
        hCenterMax(i)=max(hCenterMax(i),hi);
    end
end

% Numerical safety factors.
a=P8.selfDeflation*a0;
c=P8.crossInflation*c;
h=P8.forcingInflation*hCenterMax;

alpha=a./SM.M;
beta=zeros(8,8);
for i=1:8
    for j=1:8
        if mask(i,j)
            beta(i,j)=c(i,j)/sqrt(SM.m(i)*SM.m(j));
        end
    end
end

D=diag(alpha/2);
G=beta/2;
Gamma=D\G;
rho=max(abs(eig(Gamma)));
Mmat=D-G;

if rho<1
    v=Mmat\ones(8,1);
    w=Mmat.'\ones(8,1);
    piw=w./v;
    Pi=diag(piw);
    Q=Mmat.'*Pi+Pi*Mmat;
    minEigQ=min(eig(0.5*(Q+Q.')));
    rInf=Mmat\h;
else
    v=nan(8,1);
    w=nan(8,1);
    piw=nan(8,1);
    Pi=nan(8);
    Q=nan(8);
    minEigQ=nan;
    rInf=nan(8,1);
end

W.SM=SM;
W.points=Pts;
W.a0=a0;
W.a=a;
W.c=c;
W.alpha=alpha;
W.beta=beta;
W.D=D;
W.G=G;
W.Gamma=Gamma;
W.rho=rho;
W.M=Mmat;
W.v=v;
W.w=w;
W.pi=piw;
W.Pi=Pi;
W.Q=Q;
W.minEigQ=minEigQ;
W.hbar=h;
W.rInf=rInf;
W.maxForbidden=maxForbidden;
W.roundTripMax=roundTripMax;
W.nWitnessPoints=numel(Pts);
end
