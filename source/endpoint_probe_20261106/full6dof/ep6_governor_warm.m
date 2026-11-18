function G=ep6_governor_warm(x,u0,R,method,O)
% Same QP, weights, rows and original w2b_qp3 solver. This changes ONLY its
% numerical feasible starting-point proposal; the solver always independently
% verifies feasibility/KKT and retains its enumeration fallback.
persistent previousSupport previousMethod;
if strcmp(method,'N0'),R.DeltaN=0;end
C=ep6_affine(x,R);
if strcmp(method,'T'),ix=[1:13 16:17];else,ix=1:17;end
A=-C.G(ix,:);b=C.c(ix);Q=O.QP;Q.FeasibleStart=rg_backup_control(x,R);
% First try the updated weighted projection onto the previously supporting
% rows. This is only a proposed start, never accepted as the optimizer here.
if ~isempty(previousSupport)&&strcmp(previousMethod,method)
    W=previousSupport(:);As=bsxfun(@times,A(W,:),Q.Scale');d=b(W)-A(W,:)*u0;
    if numel(W)<=3&&rcond(As*As')>1e-12
        start=u0+Q.Scale.*(As'*((As*As')\d));
        if all(A*start<=b+1e-12*(1+abs(b)+abs(A)*abs(start))),Q.FeasibleStart=start;end
    end
end
if any(A*Q.FeasibleStart>b+1e-12*(1+abs(b)+abs(A)*abs(Q.FeasibleStart)))&&all(b>=0)
    Q.FeasibleStart=zeros(3,1);
end
[u,D]=w2b_qp3(u0,A,b,Q.Scale,Q);
if ~D.Success,error('EP6:QPNoVerifiedSolution','%s returned %s',method,D.Status);end
previousSupport=D.SupportRows;previousMethod=method;
G=struct('u',u,'C',C,'Solver',D,'phi',C.c(ix)+C.G(ix,:)*u,'rows',ix);
end
