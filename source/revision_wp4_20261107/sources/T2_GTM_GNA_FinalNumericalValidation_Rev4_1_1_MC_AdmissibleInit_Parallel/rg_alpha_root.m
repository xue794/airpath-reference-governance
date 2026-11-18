function [alpha,info] = rg_alpha_root(V,X,N,R)
%RG_ALPHA_ROOT Monotone bisection solution of N(alpha;V,X)=N.
lo=R.alphaLo; hi=R.alphaHi;
Mlo=rg_N_map(lo,V,X,R);
Mhi=rg_N_map(hi,V,X,R);
flo=Mlo.value-N;
fhi=Mhi.value-N;
info.bracket = (flo<=0) && (fhi>=0);
if ~info.bracket
    alpha=NaN;
    info.iter=0;
    info.residual=NaN;
    return;
end
for k=1:70
    mid=0.5*(lo+hi);
    Mm=rg_N_map(mid,V,X,R);
    fm=Mm.value-N;
    if fm>=0
        hi=mid;
    else
        lo=mid;
    end
end
alpha=0.5*(lo+hi);
info.iter=70;
Ma=rg_N_map(alpha,V,X,R);
info.residual=Ma.value-N;
end
