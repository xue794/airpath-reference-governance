function export_affine()
% No integration. Evaluate the original MATLAB kernel at stored or prescribed points.
here=fileparts(mfilename('fullpath'));root=fileparts(here);
addpath(fullfile(root,'revision_probe_20261106','source_snapshot','core'));
base=fullfile(root,'endpoint_probe_20261106','full6dof','V31_C60_warm_fine');
out=fullfile(here,'results'); S=load(fullfile(base,'N10.mat'));R=S.C.R;
assert(R.rhoB==.10 && R.DeltaN==10 && R.DeltaT==3);
sets={};labels={}; rates=[18 22 26];radii=[.03 .01 .003];
for nm={'N10','T'}
    v=load(fullfile(base,[nm{1} '.mat']));x=v.Xs(:,20:25);t=v.Ts(:);
    u=zeros(numel(t),3);
    for k=1:numel(t),u(k,:)=rg_backup_control(x(k,:)',R)';end
    sets{end+1}=evaluate(x,t,u,R);labels{end+1}=nm{1};
end
% Both the historical reference and the hypothetical zero-turn request use the
% same causal input before the turn (only the archived plant is available).
addpath(fileparts(base));O=S.O;O.Course=0;prefix=zeros(321,2);
for k=1:321
    t=(k-1)/20;[a,b]=ep6_reference(t,S.C,S.O);[c,d]=ep6_reference(t,S.C,O);
    prefix(k,:)=[norm(a-c,inf),norm(b-d,inf)];
end
assert(max(prefix(:))==0);
L=chol(R.PB);K=[R.kVp 0 R.kVd 0 0;0 R.kGammaP 0 R.kGammaD 0;0 0 0 0 R.kChiD];
for D=rates
    grid=(0:1100)'/20;[X,U]=guide(grid,D,R);
    sets{end+1}=evaluate(X,grid,U,R); labels{end+1}=sprintf('center%d',D);
    times=unique([0 3 5.6 14 14.65 16 16.8 18.7 20 35 40 45 50 55,3+D*[.25 .5 .75 1]]);
    [xc,uc]=guide(times,D,R);
    for radius=radii
        offs=.97*sqrt(radius)*(L\[eye(5) -eye(5)])';
        X=[];U=[];tt=[];err=[];
        for j=1:numel(times)
            o6=[offs(:,1:2),zeros(10,1),offs(:,3:5)];
            X=[X;bsxfun(@plus,xc(j,:),o6)];U=[U;bsxfun(@minus,uc(j,:),offs*K')];
            tt=[tt;repmat(times(j),10,1)];err=[err;offs];
        end
        B=evaluate(X,tt,U,R);B.offset=err;B.radius=radius;
        sets{end+1}=B;labels{end+1}=sprintf('section%d_%03d',D,round(radius*1000));
    end
end
save(fullfile(out,'affine_export.mat'),'sets','labels','R','prefix','-v7');
fprintf('OFFLINE_AFFINE_DONE %d datasets\n',numel(sets));
end

function [x,u]=guide(times,D,R)
x=zeros(numel(times),6);u=zeros(numel(times),3);
for i=1:numel(times)
    t=times(i);
    if t<=3,a=R.Vb;b=a;q=0;h=1;
    elseif t<=3+D,a=R.Vb;b=31.25;h=D;q=(t-3)/h;
    elseif t<=35,a=31.25;b=a;h=1;q=0;
    else,a=31.25;b=R.Vb;h=20;q=(t-35)/h;
    end
    x(i,1)=a+(b-a)*(10*q^3-15*q^4+6*q^5);
    x(i,4)=(b-a)*30*q^2*(1-q)^2/h;
    u(i,1)=(b-a)*60*q*(1-q)*(1-2*q)/h^2;
end
end

function B=evaluate(x,t,u,R)
n=numel(t); B.x=x;B.t=t;B.u=u;B.c=nan(n,17);B.G=nan(n,17,3);
B.phi=nan(n,17);B.h=nan(n,4);B.psi=nan(n,4);B.hRate=nan(n,3);
B.hGeom=nan(n,3);B.hForce=nan(n,2);B.hThrust=nan(n,2);B.vb=nan(n,1);
B.alpha=nan(n,1);B.slope=nan(n,1); B.T=nan(n,1); B.names={};
for k=1:n
    xf=x(k,:)';F=rg_feedforward(xf,R);
    B.hForce(k,:)=[F.hNminus F.hNplus]; z=rg_backup_state(xf,R);B.vb(k)=z'*R.PB*z;
    if ~F.rootInfo.bracket,continue;end
    a=rg_constraint_affine(xf,R);
    B.c(k,:)=a.c(1:17)';B.G(k,:,:)=a.G(1:17,:);
    B.phi(k,:)=(a.c(1:17)+a.G(1:17,:)*u(k,:)')';B.h(k,:)=a.h';B.psi(k,:)=a.psi1';
    B.hRate(k,:)=R.nuMax(:)'-abs(xf(4:6))'; B.hGeom(k,:)=[a.aGamma-R.aGammaMin,tan(R.muMax)*a.aGamma-a.aChi,tan(R.muMax)*a.aGamma+a.aChi];
    B.hThrust(k,:)=[a.T-3,57-a.T];B.alpha(k)=a.alpha;
    M=rg_N_map(a.alpha,xf(1),F.X,R);B.slope(k)=M.alpha;B.T(k)=a.T;
    B.names=a.names(1:17);
end
end
