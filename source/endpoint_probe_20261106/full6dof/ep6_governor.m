function G=ep6_governor(x,u0,R,method,O)
if strcmp(method,'N0'),R.DeltaN=0;end
C=ep6_affine(x,R);
if strcmp(method,'T'),ix=[1:13 16:17];else,ix=1:17;end
Q=O.QP;Q.FeasibleStart=rg_backup_control(x,R);
[u,D]=w2b_qp3(u0,-C.G(ix,:),C.c(ix),Q.Scale,Q);
if ~D.Success,error('EP6:QPNoVerifiedSolution','%s returned %s',method,D.Status);end
G=struct('u',u,'C',C,'Solver',D,'phi',C.c(ix)+C.G(ix,:)*u,'rows',ix);
end
