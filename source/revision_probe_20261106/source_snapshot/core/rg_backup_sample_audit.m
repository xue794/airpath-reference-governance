function A = rg_backup_sample_audit(R,Ns)
%RG_BACKUP_SAMPLE_AUDIT Dense deterministic numerical audit of K_B.
% This is a numerical module test, not a substitute for the analytical
% compact-set certificate used in the theorem.
if nargin<2, Ns=6000; end
rng(11,'twister');
[V,D]=eig(R.PB);
PinvHalf=V*diag(1./sqrt(diag(D)))*V';

minSigma=inf;
minH=inf;
minPsi=inf;
minRootMargin=inf;
worstName='';
for k=1:Ns
    d=randn(5,1); d=d/norm(d);
    rad=rand()^(1/5);
    xi=PinvHalf*(sqrt(R.rhoB)*rad*d);
    x=[R.Vb+xi(1); R.gammaB+xi(2); 0; xi(3);xi(4);xi(5)];
    C=rg_constraint_affine(x,R);
    uB=rg_backup_control(x,R);
    sig=C.c+C.G*uB;
    [ms,ii]=min(sig);
    if ms<minSigma
        minSigma=ms; worstName=C.names{ii};
    end
    minH=min(minH,min(C.h));
    minPsi=min(minPsi,min(C.psi1));
    minRootMargin=min(minRootMargin,min(C.hNminus,C.hNplus));
end
A.minSigma=minSigma; A.worstName=worstName;
A.minH=minH; A.minPsi=minPsi; A.minRootMargin=minRootMargin;
A.Ns=Ns;
end
