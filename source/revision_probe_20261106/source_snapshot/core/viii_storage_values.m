function O=viii_storage_values(X,C,P8,SM)
%VIII_STORAGE_VALUES Quadratic witness and exact F/R storage comparison.

if nargin<4 || isempty(SM), SM=viii_storage_matrices(C,P8); end
B=SM.blocks;
sc=P8.scaleVector(:);
y=X(:)./sc;

V=zeros(8,1);
for i=1:8
    yi=y(B.idx{i});
    V(i)=0.5*yi.'*SM.P{i}*yi;
end

% Exact command-filter storage at the local state.
xF=X(B.idx{6});
eF=xF(1:3); wR=xF(4:6);
nF=norm(eF);
if nF<1
    thF=asin(nF);
    psiF=1-cos(thF);
else
    psiF=inf;
end
VFexact=0.5*(wR.'*wR)+C.F.kRf*psiF+C.F.cF*(eF.'*wR);

% Exact rotational storage.
xR=X(B.idx{7});
eR=xR(1:3); eOm=xR(4:6);
nR=norm(eR);
if nR<1
    thR=asin(nR);
    psiR=1-cos(thR);
else
    psiR=inf;
end
VRexact=0.5*eOm.'*C.P.I*eOm + C.Dp.kR*psiR + ...
         C.Dp.cR*eR.'*C.P.I*eOm;

O.V=V;
O.r=sqrt(max(V,0));
O.VFexact=VFexact;
O.VRexact=VRexact;
O.VFquad=V(6);
O.VRquad=V(7);
end
