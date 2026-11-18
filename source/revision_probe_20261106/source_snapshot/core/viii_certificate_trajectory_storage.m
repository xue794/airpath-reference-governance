function T=viii_certificate_trajectory_storage(S,C,P8,W)
%VIII_CERTIFICATE_TRAJECTORY_STORAGE Root-energy history for macro or direct witness.

n=numel(S.t);

if strcmp(W.mode,'macro')
    T=viii_macro_trajectory_storage(S,C,P8,W);
    return;
end

V=zeros(n,1);
rr=zeros(n,1);

for k=1:n
    X=viii_reduced_from_full(S.t(k),S.x(k,:).',S.mode,C,P8);
    y=X./P8.scaleVector;
    V(k)=0.5*y.'*W.P*y;
    rr(k)=sqrt(max(V(k),0));
end

rdot=zeros(n,1);
rdot(1)=(rr(2)-rr(1))/(S.t(2)-S.t(1));
rdot(n)=(rr(n)-rr(n-1))/(S.t(n)-S.t(n-1));
for k=2:n-1
    rdot(k)=(rr(k+1)-rr(k-1))/(S.t(k+1)-S.t(k-1));
end

T.t=S.t;
T.V=V;
T.r=rr;
T.rdot=rdot;
end
