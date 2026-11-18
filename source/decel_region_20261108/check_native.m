function check_native()
here=fileparts(mfilename('fullpath'));root=fileparts(here);
addpath(fullfile(root,'revision_probe_20261106','source_snapshot','core'));
S=load(fullfile(root,'decel_screening_20261108','results','affine_export.mat'));
ids=find(strcmp(S.labels,'center22')|strcmp(S.labels,'section22_010'));
n=0;mx=zeros(1,4);T=[];X=[];phi=[];hB=[];
for id=ids(:)'
    B=S.sets{id};
    for i=1:numel(B.t)
        t=B.t(i);x=B.x(i,:)';[C,u,h,c]=moving_backup_affine(t,x,S.R);
        e=x([1 2 4 5 6])-c;y=C.c+C.G*u;
        expected=-e'*(S.R.AB'*S.R.PB+S.R.PB*S.R.AB)*e+S.R.kB*h;
        mx=max(mx,[norm(u-B.u(i,:)',Inf),norm(y(1:17)-B.phi(i,:)',Inf),abs(y(18)-expected),abs(h-(.01-e'*S.R.PB*e))]);
        assert(y(18)>0 && h>=0);
        T(end+1,1)=t;X(end+1,:)=x';phi(end+1,:)=y';hB(end+1,1)=h;n=n+1;
    end
end
assert(all(mx<1e-10));
save(fullfile(here,'results','native_check.mat'),'n','mx','T','X','phi','hB','-v7');
fprintf('MATLAB_NATIVE_PASS n=%d, max discrepancies %g %g %g %g\n',n,mx);
end
