function [u,D] = w2b_qp3(u0,A,b,scale,O)
%W2B_QP3 3-D strictly convex weighted projection onto A*u <= b.
%   min 0.5*sum(((u-u0)./scale).^2), all scale(i)>0.
% Base MATLAB only. Feasible / single-row exact paths, feasible-start primal
% active set, and <=3 independent-active-row enumeration fallback.
% Does NOT soften constraints, change their values, or fall back to backup.
% D.Success=false means no numerically verified KKT point was returned.
% Except a contradictory constant row, this is NOT proof of infeasibility.
% All KKT tolerances, raw residuals and normalized residuals are exposed.
if nargin<5 || isempty(O), O=w2b_options(); end
u0=double(u0(:)); b=double(b(:)); scale=double(scale(:)); A=double(A);
assert(numel(u0)==3 && numel(scale)==3 && all(isfinite(u0)) && ...
    all(isfinite(scale)) && all(scale>0),'W2B:QPInput','Invalid u0/scale.');
if isempty(A), A=zeros(0,3); end
assert(isreal(A) && isreal(b) && isreal(u0) && isreal(scale) && ...
    size(A,2)==3 && size(A,1)==numel(b) && ...
    all(isfinite(A(:))) && all(isfinite(b)), ...
    'W2B:QPInput','Expected finite real A (m-by-3) and b (m-by-1).');
assert(strcmp(O.Backend,'auto') || strcmp(O.Backend,'enumeration'), ...
    'W2B:Backend','Backend must be auto or enumeration.');
m=size(A,1); u=nan(3,1);
D=struct('Success',false,'Status','NOT_SOLVED','Path','none','Iterations',0, ...
    'Candidates',0,'RankSkipped',0,'FallbackUsed',false,'Objective',NaN, ...
    'z',nan(3,1),'Dual',zeros(m,1),'DualScaled',zeros(m,1), ...
    'ActiveRows',zeros(0,1),'SupportRows',zeros(0,1),'RowNorm',zeros(m,1), ...
    'RawResidual',nan(m,1),'RowTolerance',nan(m,1), ...
    'PrimalNormalized',Inf,'StationarityNormalized',Inf, ...
    'DualViolationNormalized',Inf,'ComplementarityNormalized',Inf, ...
    'StationarityRaw',Inf,'MinActiveSingularValue',NaN);
% z = (u-u0)./scale.  The centered objective is 0.5*z''*z.
AS=bsxfun(@times,A,scale.');
r=hypot(hypot(AS(:,1),AS(:,2)),AS(:,3));
assert(all(isfinite(r)),'W2B:Overflow','Constraint scaling overflow.');
D.RowNorm=r;
zero=(r==0);
% A true constant contradiction is mathematical infeasibility, not tolerance.
if any(zero & b<0)
    D.Status='INFEASIBLE_CONSTANT_ROW'; D.Path='constant-row'; return;
end
keep=find(~zero);
G=bsxfun(@rdivide,AS(keep,:),r(keep));
d=(b(keep)-A(keep,:)*u0)./r(keep);
assert(all(isfinite(d)),'W2B:Overflow','Centered right-hand side overflow.');
n=numel(keep); z0=zeros(3,1); mu0=zeros(n,1);
if isempty(keep)
    finish(z0,mu0,'unconstrained',0); return;
end
if strcmp(O.Backend,'auto')
    [ok,~]=kkt(z0,mu0);
    if ok, finish(z0,mu0,'feasible',0); return; end
    % Only violated single rows can have a nonzero projection multiplier.
    [~,order]=sort(d,'ascend');
    for q=1:n
        j=order(q); if d(j)>=0, continue; end
        zz=(d(j)/(G(j,:)*G(j,:).'))*G(j,:).';
        mm=zeros(n,1); mm(j)=-d(j)/(G(j,:)*G(j,:).');
        D.Candidates=D.Candidates+1;
        [ok,~]=kkt(zz,mm);
        if ok, finish(zz,mm,'single-halfspace',0); return; end
    end
    % A feasible backup is only a numerical starting point; never returned as
    % an optimizer without its own KKT certificate.
    if ~isempty(O.FeasibleStart)
        st=O.FeasibleStart(:);
        assert(numel(st)==3 && all(isfinite(st)),'W2B:Start','Invalid feasible start.');
        zz=(st-u0)./scale;
        if all(G*zz-d<=tolerances(zz))
            W=zeros(0,1);
            for it=1:O.MaxIterations
                D.Iterations=it;
                if isempty(W)
                    Q=zeros(3,0); RR=zeros(0); p=-zz;
                else
                    [Q,RR]=qr(G(W,:).',0);
                    p=-zz+Q*(Q.'*zz);
                end
                if norm(p,inf)<=O.StepTol*(1+norm(zz,inf))
                    if isempty(W), mw=zeros(0,1); else, mw=RR\(-Q.'*zz); end
                    if isempty(mw) || min(mw)>=-O.DualTol*(1+norm(zz,inf))
                        mm=zeros(n,1); mm(W)=mw;
                        [ok,~]=kkt(zz,mm);
                        if ok, finish(zz,mm,'primal-active-set',it); return; end
                        break;
                    else
                        [~,jmin]=min(mw); W(jmin)=[]; continue;
                    end
                end
                gp=G*p; slack=d-G*zz;
                inactive=true(n,1); inactive(W)=false;
                ix=find(inactive & gp>O.StepTol*(1+norm(p,inf)));
                alpha=1; blocker=0;
                if ~isempty(ix)
                    aa=max(0,slack(ix))./gp(ix);
                    [aval,jmin]=min(aa);
                    if aval<1, alpha=aval; blocker=ix(jmin); end
                end
                zz=zz+alpha*p;
                if any(G*zz-d>10*tolerances(zz)), break; end
                if blocker>0
                    WW=[W;blocker]; sv=svd(G(WW,:));
                    if numel(WW)>3 || min(sv)<=O.RankTol*max(sv)
                        % Never regularize/drop a hard row to force progress.
                        break;
                    end
                    W=WW;
                end
            end
        end
    end
    D.FallbackUsed=true;
end
% Reference/full fallback. In R^3 a polyhedral normal cone admits a KKT
% support using at most three linearly independent row normals. Degenerate
% binding rows may have zero multipliers and are still checked for feasibility.
best=[]; bestmu=[]; bestJ=Inf;
[ok,~]=kkt(z0,mu0);
if ok
    finish(z0,mu0,'enumeration-zero',D.Iterations); return;
end
for k=1:min(3,n)
    sets=nchoosek(1:n,k);
    for q=1:size(sets,1)
        W=sets(q,:).'; B=G(W,:);
        sv=svd(B);
        if min(sv)<=O.RankTol*max(sv)
            D.RankSkipped=D.RankSkipped+1; continue;
        end
        % QR avoids squaring the active matrix condition number.
        [Q,RR]=qr(B.',0);
        yy=RR.'\d(W); zz=Q*yy;
        mw=RR\(-yy); mm=zeros(n,1); mm(W)=mw;
        D.Candidates=D.Candidates+1;
        [ok,~]=kkt(zz,mm);
        if ok
            jj=0.5*(zz.'*zz);
            if ~O.Exhaustive
                finish(zz,mm,'enumeration',D.Iterations); return;
            end
            if jj<bestJ, best=zz; bestmu=mm; bestJ=jj; end
        end
    end
end
if ~isempty(best)
    finish(best,bestmu,'enumeration-exhaustive',D.Iterations); return;
end
D.Status='NO_VERIFIED_KKT_POINT'; D.Path='enumeration-unresolved';

    function tt=tolerances(zz)
        tt=O.FeasTol+O.RelTol*(1+abs(d)+abs(G)*abs(zz));
    end
    function [pass,K]=kkt(zz,mm)
        K=struct('primal',Inf,'stationarity',Inf,'dual',Inf,'complementarity',Inf);
        if any(~isfinite(zz)) || any(~isfinite(mm)), pass=false; return; end
        res=G*zz-d; st=zz+G.'*mm;
        K.primal=max([0;res./(1+abs(d)+abs(G)*abs(zz))]);
        K.stationarity=norm(st,inf)/(1+norm(zz,inf)+norm(G.'*mm,inf));
        K.dual=max([0;-mm])/(1+norm(mm,inf));
        K.complementarity=norm(mm.*res,inf)/(1+norm(mm,inf)*(1+norm(d,inf))+zz.'*zz);
        pass=all(res<=tolerances(zz)) && ...
            all(mm>=-O.DualTol*(1+norm(zz,inf))) && ...
            K.stationarity<=O.KKTTol && K.complementarity<=O.KKTTol;
    end
    function finish(zz,mm,pathName,it)
        [ok,K]=kkt(zz,mm);
        uu=u0+scale.*zz;
        orig=A*uu-b;
        rowTol=zeros(m,1);
        rowTol(keep)=r(keep).*tolerances(zz)+ ...
            64*eps*(1+abs(b(keep))+abs(A(keep,:))*abs(uu));
        rowTol(zero)=0;
        D.RawResidual=orig; D.RowTolerance=rowTol;
        if ~ok || any(orig>rowTol)
            D.Status='ORIGINAL_ROW_CHECK_FAILED'; D.Path=pathName; return;
        end
        u=uu; D.Success=true; D.Status='OPTIMAL_KKT_VERIFIED';
        D.Path=pathName; D.Iterations=it; D.Objective=0.5*(zz.'*zz);
        D.z=zz; D.DualScaled(keep)=mm; D.Dual(keep)=mm./r(keep);
        D.ActiveRows=keep(abs(G*zz-d)<=10*tolerances(zz));
        D.SupportRows=keep(mm>10*O.DualTol*(1+norm(mm,inf)));
        D.PrimalNormalized=K.primal; D.StationarityNormalized=K.stationarity;
        D.DualViolationNormalized=K.dual; D.ComplementarityNormalized=K.complementarity;
        D.StationarityRaw=norm((u-u0)./(scale.^2)+A.'*D.Dual,inf);
        if ~isempty(D.SupportRows)
            B=bsxfun(@rdivide,AS(D.SupportRows,:),r(D.SupportRows));
            sv=svd(B); D.MinActiveSingularValue=min(sv);
        end
    end
end
