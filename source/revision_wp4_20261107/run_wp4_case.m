function run_wp4_case(mode)
% Exactly one frozen missing cell. Original RHS, solver and tolerances.
H=fileparts(mfilename('fullpath'));addpath(fullfile(H,'sources','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel'),'-begin');
addpath(fullfile(H,'sources','transparent'),'-begin');
assert(any(strcmp(mode,{'hold','noisy'})),'WP4:Mode','Only the two missing cells are authorized here.');
pre=jsondecode(fileread(fullfile(H,'checks','preflight.json')));assert(pre.Passed);
key=['Transparent_' mode];target=fullfile(H,'runs',[key '.mat']);assert(exist(target,'file')~=2,'WP4:AlreadyRun','Refusing to overwrite an existing trajectory.');
F=load(fullfile(H,'inputs',[key '_frozen.mat']));I=F.I;C=I.C;K=I.K;x0=I.x0;
assert(C.R.rhoB==.10 && C.RelTol==1e-8 && C.AbsTol==1e-10 && C.MaxStep==.01 && C.sample_dt==.05);
tgrid=(0:C.sample_dt:K.T).';savedT=0;savedX=x0.';lastReport=-10;attemptT=0;rhsCalls=0;wall=tic;
status='running';errorId='';errorText='';lastwarn('');
fprintf('WP4_RUN_START %s T=%.2f dt=%.3f MaxStep=%.3f no retuning\n',key,K.T,C.sample_dt,C.MaxStep);
opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep,'OutputFcn',@observe);
try
 [t,x]=ode45(@rhs,tgrid,x0,opts);
 if abs(t(end)-K.T)<1e-10,status='completed';else,status='integrator_stopped_early';end
catch ME
 t=savedT;x=savedX;status=classify(ME.identifier);errorId=ME.identifier;errorText=getReport(ME,'extended');
 fprintf(2,'WP4_RUN_EXCEPTION %s last returned t %.12g attempted t %.12g %s\n',key,t(end),attemptT,errorId);
end
elapsed=toc(wall);[warningText,warningId]=lastwarn;
Baseline=struct('t',t,'x',x,'mode',mode,'Case',K,'I0',I.I0,'C',C,'Method','Transparent', ...
 'RunStatus',status,'ExpectedDuration',K.T,'IntegrationWall_s',elapsed,'RHSCalls',rhsCalls,'ErrorId',errorId,'ErrorText',errorText);
save(target,'Baseline','-v7'); % Preserve the raw trajectory before any output reconstruction.
try
 Baseline.Y=cl6_outputs_transparent_minimal_R2018a(t,x,mode,C);
 Baseline.ablation=struct('name','transparent_governor','certificationEnabled',false,'lambdaForced',1,'retuned',false,'mode',mode);
 save(target,'Baseline','-v7');
catch ME
 errorId=ME.identifier;errorText=getReport(ME,'extended');status='postprocessing_error_raw_trajectory_retained';
end
R=struct('Case',mode,'Method','Transparent','Status',status,'Completed',strcmp(status,'completed'),'ExpectedDuration_s',K.T, ...
 'LastReturnedTime_s',t(end),'Samples',numel(t),'ErrorId',errorId,'ErrorText',errorText,'LastAttemptedRHSTime_s',attemptT, ...
 'IntegrationWall_s',elapsed,'RHSCalls',rhsCalls,'WarningId',warningId,'WarningText',warningText, ...
 'NewIntegrationUnits',1,'NoRetuning',true,'File',target,'MATLAB',version);
f=fopen(fullfile(H,'runs',[key '_status.json']),'w');fprintf(f,'%s',jsonencode(R));fclose(f);
fprintf('WP4_RUN_FINISHED %s %s n=%d t=%.6f wall=%.1fs\n',key,status,numel(t),t(end),elapsed);

 function dx=rhs(tt,xx)
  attemptT=tt;rhsCalls=rhsCalls+1;
  dx=cl6_rhs_transparent_R2018a(tt,xx,mode,C);
 end
 function stop=observe(tt,xx,flag)
  stop=0;
  if isempty(flag)
   for ii=1:numel(tt)
    if tt(ii)>savedT(end)+1e-12,savedT(end+1,1)=tt(ii);savedX(end+1,:)=xx(:,ii).';end %#ok<AGROW>
   end
   if ~isempty(tt) && tt(end)-lastReport>=10-1e-8
    lastReport=tt(end);progress=struct('case',mode,'last_saved_s',tt(end),'wall_s',toc(wall),'rhs_calls',rhsCalls,'state','running');
    f=fopen(fullfile(H,'runs',[key '_progress.json']),'w');fprintf(f,'%s',jsonencode(progress));fclose(f);
    fprintf('WP4_PROGRESS %s t=%.2f/%.2f wall=%.1fs\n',key,tt(end),K.T,toc(wall));
    save(fullfile(H,'runs',[key '_checkpoint.mat']),'savedT','savedX','-v7');
   end
  end
 end
end

function label=classify(id)
switch id
 case 'rg_raw_command_chart:singular',label='raw_command_chart_exit';
 case 'rg_constraint_affine:root',label='reference_inversion_domain_exit';
 case 'viA_feedback_command:root',label='feedback_inversion_domain_exit';
 case 'transparent_governor:kappa',label='thrust_projection_domain_exit';
 otherwise,label='unclassified_exception';
end
end
