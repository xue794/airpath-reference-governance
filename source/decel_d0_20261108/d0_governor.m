function G=d0_governor(t,x,u0,R,method,O)
persistent previousSupport previousMethod
if any(strcmp(method,{'F','D','A'}))
    [C,ub]=moving_backup_affine(t,x,R);
else
    C=cp6_affine(x,R);ub=rg_backup_control(x,R);
end
switch method
 case {'O','F','A'}, ix=1:18;
 case 'D', ix=[1:13,16:18];
 case {'N10','Am'}, ix=1:17;
 case 'T', ix=[1:13,16:17];
 otherwise, error('Unexpected D0 arm');
end
cc=C.c(ix); GG=C.G(ix,:);
if any(strcmp(method,{'A','Am'}))
    Ac=cp_aoa_rows(x,R,'Am',C.alpha);
    cc(14:15)=Ac.c;GG(14:15,:)=Ac.G;
end
A=-GG;b=cc;Q=O.QP;Q.FeasibleStart=ub;
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
G=struct('u',u,'C',C,'Solver',D,'phi',cc+GG*u,'rows',ix,'cc',cc,'GG',GG,'backup',ub);
end
