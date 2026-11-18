function G=cp6_governor_warm(x,u0,R,method,O)
% Same QP, weights, rows and original w2b_qp3 solver. This changes ONLY its
% numerical feasible starting-point proposal; the solver always independently
% verifies feasibility/KKT and retains its enumeration fallback.
persistent previousSupport previousMethod;
if strcmp(method,'N0')||strcmp(method,'A0'),R.DeltaN=0;end
C=cp6_affine(x,R);
if strcmp(method,'A0')||strcmp(method,'Am')
    ix=1:17;Ac=cp_aoa_rows(x,R,method,C.alpha);
    cc=C.c(ix);GG=C.G(ix,:);cc(14:15)=Ac.c;GG(14:15,:)=Ac.G;
elseif strcmp(method,'T')
    ix=[1:13 16:17];cc=C.c(ix);GG=C.G(ix,:);
else
    ix=1:17;cc=C.c(ix);GG=C.G(ix,:);
end
A=-GG;b=cc;Q=O.QP;Q.FeasibleStart=rg_backup_control(x,R);
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
if ~D.Success,error('CP6:QPNoVerifiedSolution','%s returned %s',method,D.Status);end
previousSupport=D.SupportRows;previousMethod=method;
G=struct('u',u,'C',C,'Solver',D,'phi',cc+GG*u,'rows',ix);
end
