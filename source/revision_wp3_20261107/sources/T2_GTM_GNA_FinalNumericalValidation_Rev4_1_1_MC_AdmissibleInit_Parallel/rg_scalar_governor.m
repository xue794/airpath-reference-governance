function O = rg_scalar_governor(x,u0,R)
%RG_SCALAR_GOVERNOR Certified interpolation u*=uB+lambda*(u0-uB).
C=rg_constraint_affine(x,R);
uB=rg_backup_control(x,R);
sigma=C.c+C.G*uB;
if min(sigma)<-1e-9
    error('rg_scalar_governor:backup','Backup command lost feasibility: min sigma=%g',min(sigma));
end
sigma=max(sigma,0);

d=u0-uB;
eta=-(C.G*d); % a_j'd for phi=b-a'u

lambdaJ=ones(size(sigma));
idx=eta>sigma;
lambdaJ(idx)=sigma(idx)./eta(idx);
lambda=min(lambdaJ);
lambda=max(0,min(1,lambda));

u=uB+lambda*d;
phi=C.c+C.G*u;

O.u=u; O.lambda=lambda; O.lambdaJ=lambdaJ;
O.u0=u0; O.uB=uB; O.sigma=sigma; O.phi=phi;
O.C=C;
end
