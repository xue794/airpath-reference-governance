function run_turn_crosscheck()
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
addpath(fullfile(root,'revision_probe_20261106','source_snapshot','core'));
R=rg_params();R.Vb=34;R.rhoB=.005;turn=.15;
half=sqrt(R.rhoB*diag(inv(R.PB)));
Z=[zeros(1,5);diag(.4*half);-diag(.4*half)];Rows=cell(size(Z,1),1);
for k=1:size(Z,1)
 z=Z(k,:).';x=[R.Vb+z(1);z(2);0;z(3);z(4);turn+z(5)];
 A=rg_constraint_affine(x,R);u=rg_backup_control(x,R);u(3)=-R.kChiD*z(5);
 Rows{k}=struct('z',z,'x',x,'u',u,'phi',A.c(1:17)+A.G(1:17,:)*u, ...
   'Tff',A.T,'alpha',A.alpha,'hNminus',A.hNminus,'hNplus',A.hNplus);
end
fid=fopen(fullfile(here,'turn_matlab_crosscheck.json'),'w');fprintf(fid,'%s',jsonencode(Rows));fclose(fid);
fprintf('Saved %d relative-turn crosscheck points\n',numel(Rows));
end
