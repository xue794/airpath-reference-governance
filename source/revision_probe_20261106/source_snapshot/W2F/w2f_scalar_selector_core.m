function [u,D] = w2f_scalar_selector_core(u0,C,uB,scale)
%W2F_SCALAR_SELECTOR_CORE Frozen scalar selection with prebuilt C and uB.
% Used only to isolate selector cost. Algebra matches rg_scalar_governor.
u0=u0(:);uB=uB(:);scale=scale(:);
sigma=C.c+C.G*uB;
if min(sigma)<-1e-9
    error('W2F:Backup','Backup command infeasible in selector benchmark: %.17g',min(sigma));
end
sigma=max(sigma,0);
d=u0-uB;
eta=-(C.G*d);
lambdaJ=ones(size(sigma));
idx=eta>sigma;
lambdaJ(idx)=sigma(idx)./eta(idx);
lambda=min(lambdaJ);
lambda=max(0,min(1,lambda));
u=uB+lambda*d;
phi=C.c+C.G*u;
D=struct('u',u,'lambda',lambda,'lambdaJ',lambdaJ,'sigma',sigma,'phi',phi, ...
    'Objective',0.5*sum(((u-u0)./scale).^2), ...
    'CorrectionNorm',norm((u-u0)./scale));
end
