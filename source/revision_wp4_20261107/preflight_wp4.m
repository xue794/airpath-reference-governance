function preflight_wp4()
% Algebraic identity checks and full input freeze only. NO integrations.
H=fileparts(mfilename('fullpath'));core=fullfile(H,'sources','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel');
addpath(core,'-begin');addpath(fullfile(H,'sources','P1_Aerospace_W2C_R2018a_v1_0'),'-begin');addpath(fullfile(H,'sources','transparent'),'-begin');
C=cl6_params();C.R.rhoB=.10;O=w2c_options();
Q=load(fullfile(fileparts(H),'revision_wp3_20261107','results','Proposed_hold_export.mat'),'C');assert(isequaln(C,Q.C));
required={'cl6_params','cl6_case_info','cl6_initial_state','cl6_reference','cl6_wind_profile','wind_estimator_test_noise', ...
 'cl6_rhs_transparent_R2018a','cl6_control_snapshot_transparent_R2018a','rg_constraint_affine','rg_alpha_root','viA_feedback_command','t2_gna_6dof_quat','vii_actuator_derivative'};
paths=cell(size(required));for i=1:numel(required),paths{i}=which(required{i});assert(strncmp(paths{i},fullfile(H,'sources'),numel(fullfile(H,'sources'))));end
cases={'hold','noisy'};Tests={};
for j=1:2
 mode=cases{j};K=cl6_case_info(mode,C);[x0,I0]=cl6_initial_state(mode,C);
 old=load(fullfile(fileparts(H),'revision_wp3_20261107','results',['Proposed_' mode '_export.mat']));
 assert(max(abs(x0-old.x(1,:).'))<1e-12);
 idx=unique(round(linspace(1,numel(old.t),12)));err=zeros(numel(idx),3);holdFeasible=true;
 for k=1:numel(idx)
  ix=idx(k);t=old.t(ix);x=old.x(ix,:).';A=cl6_control_snapshot_transparent_R2018a(t,x,mode,C);B=w2c_control_snapshot(t,x,mode,C,'Proposed',O);
  fA=cl6_rhs_transparent_R2018a(t,x,mode,C);fB=w2c_rhs(t,x,mode,C,'Proposed',O);
  keep=[1:22,26:36];err(k,:)=[max(abs(fA(keep)-fB(keep))),max(abs((fA(23:25)-fB(23:25))-(A.Gov.u-B.Gov.u))),max(abs(A.Gov.u-B.u0))];
  assert(max(err(k,:))<1e-12);
  if strcmp(mode,'hold'),holdFeasible=holdFeasible && max(abs(fA-fB))<1e-12;end
 end
 assert(holdFeasible);
 src={'cl6_reference.m','cl6_wind_profile.m','wind_estimator_test_wind.m','wind_estimator_test_noise.m','wind_estimator_params.m'};
 signalSource=struct();for k=1:numel(src),[~,nm]=fileparts(src{k});signalSource.(nm)=fileread(fullfile(core,src{k}));end
 I=struct('Case',mode,'Method','Transparent','C',C,'K',K,'x0',x0,'I0',I0,'MATLAB',version,'Paths',{paths}, ...
  'RequiredFunctions',{required},'Signals',signalSource,'SelectedIndices',idx,'IdentityErrors',err,'HoldRHSMatches',holdFeasible,'Integrations',0);
 save(fullfile(H,'inputs',['Transparent_' mode '_frozen.mat']),'I','-v7');
 fid=fopen(fullfile(H,'inputs',['Transparent_' mode '_frozen.json']),'w');fprintf(fid,'%s',jsonencode(I));fclose(fid);
 Tests{j}=struct('Case',mode,'Passed',true,'InitMatchesOriginal',true,'MaxDownstreamDifference',max(err(:,1)), ...
  'MaxDerivativeSelectorIdentityError',max(err(:,2)),'MaxAppliedU0Difference',max(err(:,3)));
end
Report=struct('Tests',[Tests{:}],'ParametersMatchWP3',true,'Integrations',0,'MATLAB',version,'Passed',true);
fid=fopen(fullfile(H,'checks','preflight.json'),'w');fprintf(fid,'%s',jsonencode(Report));fclose(fid);
fprintf('WP4_PREFLIGHT_PASS no integrations\n');
end
