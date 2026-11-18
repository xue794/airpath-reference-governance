function T=viii_trajectory_storage(S,C,P8,W)
%VIII_TRAJECTORY_STORAGE Evaluate reduced subsystem/composite storage time histories.

n=numel(S.t);
V=zeros(n,8);
r=zeros(n,8);
Xn=zeros(n,8);

for k=1:n
    [X,~]=viii_reduced_from_full(S.t(k),S.x(k,:).',S.mode,C,P8);
    O=viii_storage_values(X,C,P8,W.SM);
    V(k,:)=O.V.';
    r(k,:)=O.r.';

    B=W.SM.blocks;
    y=X./P8.scaleVector;
    for i=1:8
        Xn(k,i)=norm(y(B.idx{i}));
    end
end

Vcomp=V*W.pi;
T.V=V;
T.r=r;
T.Xnorm=Xn;
T.Vcomp=Vcomp;
end
