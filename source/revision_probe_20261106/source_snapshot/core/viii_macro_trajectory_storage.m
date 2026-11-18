function T=viii_macro_trajectory_storage(S,C,P8,W)
%VIII_MACRO_TRAJECTORY_STORAGE Macro root-energy histories on a real trajectory.

n=numel(S.t);
r=zeros(n,2);
V=zeros(n,2);

for k=1:n
    X=viii_reduced_from_full(S.t(k),S.x(k,:).',S.mode,C,P8);
    y=X./P8.scaleVector;

    y1=y(W.idx1);
    y2=y(W.idx2);

    V(k,1)=0.5*y1.'*W.G1.P*y1;
    V(k,2)=0.5*y2.'*W.G2.P*y2;

    r(k,:)=sqrt(max(V(k,:),0));
end

% Central finite-difference root derivatives on the stored uniform output grid.
rdot=zeros(n,2);
for j=1:2
    rdot(1,j)=(r(2,j)-r(1,j))/(S.t(2)-S.t(1));
    rdot(n,j)=(r(n,j)-r(n-1,j))/(S.t(n)-S.t(n-1));
    for k=2:n-1
        rdot(k,j)=(r(k+1,j)-r(k-1,j))/(S.t(k+1)-S.t(k-1));
    end
end

T.t=S.t;
T.V=V;
T.r=r;
T.rdot=rdot;
end
