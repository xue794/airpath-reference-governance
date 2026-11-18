function G=w2c_scalar_from_constraints(x,u0,R,C,O)
%W2C_SCALAR_FROM_CONSTRAINTS Exact scalar selection given preconstructed C; timing/development only.
uB=rg_backup_control(x,R);sigma=C.c+C.G*uB;if min(sigma)<-1e-9,error('W2C:Backup','Backup infeasible.');end
sigma=max(sigma,0);d=u0-uB;eta=-(C.G*d);lambdaJ=ones(size(sigma));idx=eta>sigma;lambdaJ(idx)=sigma(idx)./eta(idx);
lambda=min(lambdaJ);lambda=max(0,min(1,lambda));u=uB+lambda*d;phi=C.c+C.G*u;
G=struct('u',u,'uB',uB,'lambda',lambda,'lambdaJ',lambdaJ,'phi',phi, ...
 'Objective',0.5*sum(((u-u0(:))./O.QP.Scale).^2));
end
