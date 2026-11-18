function run_backup_probe()
% Read-only audit of frozen constraints on retained round-two trajectories.
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
addpath(fullfile(root,'revision_probe_20261106','source_snapshot','core'));
R=rg_params();R.rhoB=.1;centers=[31 32 34 38.1];
Report.scope='Saved-state audit only; violations disprove pointwise backup feasibility, positive samples do not certify a region.';
Report.centers=centers;Report.originalR=R;
Report.cases={};
for v=[31 32]
 p=fullfile(root,'endpoint_probe_20261106','full6dof',sprintf('V%d_C60_warm_fine',v),'N10.mat');
 S=load(p); X=S.Xs(:,20:25);t=S.Ts;n=numel(t);Phi=zeros(n,17,numel(centers));VB=zeros(n,numel(centers));U=zeros(n,3,numel(centers));
 for k=1:n
  a=rg_constraint_affine(X(k,:).',R);
  for j=1:numel(centers)
   Rj=R;Rj.Vb=centers(j);z=rg_backup_state(X(k,:).',Rj);VB(k,j)=z.'*Rj.PB*z;
   u=rg_backup_control(X(k,:).',Rj);U(k,:,j)=u.';Phi(k,:,j)=(a.c(1:17)+a.G(1:17,:)*u).';
  end
 end
 B=struct('speed_request',v,'samples',n,'source',p); B.centers=cell(1,numel(centers));
 for j=1:numel(centers)
  [mn,idx]=min(Phi(:,:,j),[],1);first=find(VB(:,j)>.1+1e-12,1);
  q=struct('Vb',centers(j),'covered_samples_rho01',sum(VB(:,j)<=.1+1e-12),'minVB',min(VB(:,j)),'maxVB',max(VB(:,j)), ...
   'required_rho_to_cover_saved_trajectory',max(VB(:,j)),'first_outside_time_s',NaN,'backup_all17_pass_samples',sum(all(Phi(:,:,j)>=-1e-9,2)), ...
   'row_min',mn,'row_negative_count',sum(Phi(:,:,j)<-1e-9,1),'row_min_time_s',reshape(t(idx),1,[]),'row_min_x',X(idx,:), 'row_min_u',U(idx,:,j));
   if ~isempty(first),q.first_outside_time_s=t(first);end
   B.centers{j}=q;
 end
 Report.cases{end+1}=B;
 save(fullfile(here,sprintf('V%d_backup_saved_audit.mat',v)),'X','t','Phi','VB','U','R','centers');
  fprintf('V%d original KB coverage %d/%d; max VB %.6g; backup all17 pass %d/%d\n',v,B.centers{4}.covered_samples_rho01,n,B.centers{4}.maxVB,B.centers{4}.backup_all17_pass_samples,n);
end
Report.row_names=a.names(1:17);
% Equilibrium centers and the disclosed same-state endpoint witness.
Report.equilibria=cell(1,numel(centers));
for j=1:numel(centers)
 Rj=R;Rj.Vb=centers(j);x=[centers(j);0;0;0;0;0];a=rg_constraint_affine(x,Rj);
 Report.equilibria{j}=struct('Vb',centers(j),'x',x,'Tff',a.T,'alpha_deg',rad2deg(a.alpha),'hNminus',a.hNminus,'hNplus',a.hNplus,'phi',a.c(1:17));
end
save(fullfile(here,'backup_probe.mat'),'Report');fid=fopen(fullfile(here,'backup_probe.json'),'w');fprintf(fid,'%s',jsonencode(Report));fclose(fid);
end
