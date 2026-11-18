function Report=run_p1_aerospace_w2c(sourceRoot,varargin)
%RUN_P1_AEROSPACE_W2C One-click deterministic full-6DOF Proposed vs Matched-QP comparison.
% MATLAB R2018a. Frozen source is read-only. No Monte Carlo is run here.
p=inputParser;p.addRequired('sourceRoot',@ischar);p.addParameter('UseParallel',true,@(x)islogical(x)||ismember(x,[0 1]));p.addParameter('NumWorkers',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);p.addParameter('SaveFullTrajectories',true,@(x)islogical(x)||ismember(x,[0 1]));p.parse(sourceRoot,varargin{:});A=p.Results;
packageRoot=fileparts(mfilename('fullpath'));oldPath=path;oldPwd=pwd;cd(packageRoot);addpath(sourceRoot,'-begin');addpath(packageRoot,'-begin');rehash path;
cleanupPath=onCleanup(@()local_restore(oldPath,oldPwd)); %#ok<NASGU>
sourceRoot=local_normalize_source(sourceRoot);O=w2c_options();O.UseParallel=logical(A.UseParallel);O.RequestedWorkers=round(A.NumWorkers);O.SaveFullTrajectories=logical(A.SaveFullTrajectories);
stamp=datestr(now,'yyyymmdd_HHMMSS');outdir=[sourceRoot '_W2C_' stamp];if exist(outdir,'dir')==7,error('W2C:OutputExists','Output already exists: %s',outdir);end;mkdir(outdir);
logfile=fullfile(outdir,'W2C_console_log.txt');diary(logfile);cleanupDiary=onCleanup(@()local_diary_off()); %#ok<NASGU>
fprintf('\n============================================================\nPaper 1 -> Aerospace | W2-C deterministic full-6DOF comparison\n============================================================\n');
fprintf('Version : %s\nMATLAB  : %s\nSource  : %s\nPackage : %s\nOutput  : %s\n',O.Version,version,sourceRoot,packageRoot,outdir);
fprintf('Methods : Proposed scalar interpolation vs backup-set-matched CBF/HOCBF-QP\nModes   : hold, varying, noisy\nrhoB    : %.2f (in memory only)\n',O.RhoB);
fprintf('No frozen source file will be modified. No Monte Carlo is run.\n');
Preflight=w2c_preflight(sourceRoot,outdir,O);fprintf('PREFLIGHT PASS: exact frozen source, 18 hard rows, locked QP metric.\n');
methods=O.Methods;modes=O.Modes;N=numel(methods)*numel(modes);jobMode=cell(N,1);jobMethod=cell(N,1);q=0;for im=1:numel(modes),for it=1:numel(methods),q=q+1;jobMode{q}=modes{im};jobMethod{q}=methods{it};end;end
fprintf('\nPART A -- six paired deterministic trajectories\n');for j=1:N,fprintf('  job %d: %-9s | %s\n',j,jobMethod{j},jobMode{j});end
[P,pool]=w2c_parallel_prepare(O.UseParallel,O.RequestedWorkers,N,sourceRoot,packageRoot,outdir);
Jobs=cell(N,1);
try
 if P.Enabled
  parfor j=1:N,Jobs{j}=w2c_run_job(jobMode{j},jobMethod{j},O);end
 else
  for j=1:N,Jobs{j}=w2c_run_job(jobMode{j},jobMethod{j},O);end
 end
catch ME
 w2c_parallel_finish(pool,P,sourceRoot,packageRoot,outdir);pool=[];rethrow(ME);
end
PoolFinish=w2c_parallel_finish(pool,P,sourceRoot,packageRoot,outdir);pool=[];
fprintf('%s\n',PoolFinish.Message);
Cases=struct();Entries=cell(N,1);allOK=true;
for j=1:N
 J=Jobs{j};Entries{j}=J.Entry;key=[J.Method '_' J.Mode];key=strrep(key,'-','_');
 if J.Success
  Cases.(key)=J.S;E=J.Entry;fprintf('[%-9s | %-7s] COMPLETED: RMS=%.6f m, peak=%.6f m, intervention=%.2f%%, maxVB/rho=%.9f\n',J.Method,J.Mode,E.PosRMS_Full_m,E.PosPeak_Full_m,E.GovInterventionPct,E.MaxVB_over_rhoB);
  if O.SaveFullTrajectories,S=J.S;save(fullfile(outdir,[key '.mat']),'S','-v7.3');end %#ok<NASGU>
 else
  allOK=false;fprintf(2,'[%-9s | %-7s] FAILED: %s -- %s\n',J.Method,J.Mode,J.ErrorId,J.ErrorMessage);local_write_error(fullfile(outdir,[key '_ERROR.txt']),J);
 end
end
% Export summary with union of fields (failed rows are expanded to the complete schema when possible).
Template=[];for j=1:N,if Jobs{j}.Success,Template=Jobs{j}.Entry;break;end;end
for j=1:N,if ~Jobs{j}.Success && ~isempty(Template),Entries{j}=local_expand_failed(Jobs{j}.Entry,Template);end;end
SummaryTable=struct2table([Entries{:}]);writetable(SummaryTable,fullfile(outdir,'W2C_case_summary.csv'));
PairRows=cell(0,1);RoleRows=cell(0,1);PhiRows=cell(0,1);MarginRows=cell(0,1);
if allOK
 fprintf('\nPART B -- paired own-trajectory audit\n');
 for im=1:numel(modes)
  SP=Cases.(['Proposed_' modes{im}]);SQ=Cases.(['MatchedQP_' modes{im}]);EP=w2c_summarize_case(SP,O);EQ=w2c_summarize_case(SQ,O);pr=w2c_pairwise_row(SP,SQ,EP,EQ);PairRows{end+1}=pr; %#ok<AGROW>
  fprintf('%-7s: QP-Proposed dRMS=%+.6f m, dPeak=%+.6f m, QP thrust-support samples=%d, nonbackup-physical support=%d\n',modes{im},pr.QP_minus_Proposed_RMS_m,pr.QP_minus_Proposed_Peak_m,pr.QPThrustSupportSamples,pr.QPNonbackupPhysicalSupportSamples);
 end
 for j=1:N
  S=Cases.([Jobs{j}.Method '_' Jobs{j}.Mode]);rr=w2c_role_rows(S,O);RoleRows=[RoleRows;rr(:)]; %#ok<AGROW>
  [pp,mm]=w2c_margin_rows(S);PhiRows=[PhiRows;pp(:)];MarginRows=[MarginRows;mm(:)]; %#ok<AGROW>
 end
 writetable(struct2table([PairRows{:}]),fullfile(outdir,'W2C_pairwise_comparison.csv'));
 writetable(struct2table([RoleRows{:}]),fullfile(outdir,'W2C_constraint_roles.csv'));
 writetable(struct2table([PhiRows{:}]),fullfile(outdir,'W2C_affine_row_minima.csv'));
 writetable(struct2table([MarginRows{:}]),fullfile(outdir,'W2C_reference_margin_minima.csv'));
else
 fprintf(2,'PART B skipped because at least one trajectory failed.\n');
end
% Close parallel pool before timing to keep development timing isolated/client-serial.
fprintf('\nPART C -- development-only client-serial timing probe\n');Timing=w2c_timing_probe(Cases,O);
if strcmp(Timing.Status,'COMPLETED')
 fprintf('Timing (NOT publication benchmark): scalar selector median=%.6f ms, QP selector=%.6f ms, ratio=%.3f; whole governor ratio=%.3f\n',Timing.ScalarSelectionMedian_ms,Timing.QPSelectionMedian_ms,Timing.QP_over_ScalarSelection,Timing.QP_over_ScalarWhole);
else,fprintf('Timing skipped.\n');end
writetable(struct2table(Timing),fullfile(outdir,'W2C_timing_development.csv'));
% Gates.
CorePass=allOK;if allOK
 for j=1:N,E=Jobs{j}.Entry;CorePass=CorePass && E.ReferenceGatePass && E.ActuatorBoxPass && E.QPKKTGatePass;end
 pairModes=cellfun(@(x)x.Mode,PairRows,'UniformOutput',false);ih=find(strcmp(pairModes,'hold'),1);assert(~isempty(ih),'W2C:HoldPair','Hold pair missing.');holdPair=PairRows{ih};CorePass=CorePass && holdPair.InitialStateInfDiff<=eps && holdPair.MaxStateGridInfDiff<=O.HoldStateMatchTol;
 stressSupport=0;for k=1:numel(PairRows),if strcmp(PairRows{k}.Mode,'varying')||strcmp(PairRows{k}.Mode,'noisy'),stressSupport=stressSupport+PairRows{k}.QPNonbackupPhysicalSupportSamples;end;end
 if stressSupport>0,EvidenceStatus='EXISTING_STRESS_CASE_HAS_QP_NONBACKUP_PHYSICAL_SUPPORT';else,EvidenceStatus='TARGETED_CASE_STILL_NEEDED_FOR_COUPLED_CONSTRAINT_EVIDENCE';end
else,EvidenceStatus='NOT_EVALUATED_DUE_TO_TRAJECTORY_FAILURE';end
if CorePass,Decision='PASS_W2C_CORE_DETERMINISTIC';else,Decision='HOLD_W2C_CORE_FAILURE';end
fprintf('\n============================================================\nW2-C DECISION: %s\nEVIDENCE STATUS: %s\n============================================================\n',Decision,EvidenceStatus);
fprintf('This is a deterministic closed-loop gate; it is NOT Monte Carlo evidence and NOT a superiority verdict.\n');
Report=struct('Timestamp',datestr(now,31),'Version',O.Version,'SourceRoot',sourceRoot,'PackageRoot',packageRoot,'OutputDir',outdir,'Options',O,'Preflight',Preflight,'Parallel',P,'ParallelFinish',PoolFinish,'Jobs',{Jobs},'Cases',Cases,'PairRows',{PairRows},'RoleRows',{RoleRows},'PhiRows',{PhiRows},'MarginRows',{MarginRows},'Timing',Timing,'Decision',Decision,'EvidenceStatus',EvidenceStatus);
save(fullfile(outdir,'W2C_Report.mat'),'Report','-v7.3');local_write_readme(fullfile(outdir,'W2C_README_RESULTS.txt'),Report);
diary off;
zipPath=[outdir '.zip'];try,zip(zipPath,{'*'},outdir);fprintf('Output ZIP: %s\n',zipPath);catch ME,fprintf(2,'ZIP creation failed: %s\n',ME.message);end
fprintf('Output directory: %s\n',outdir);
end
function p=local_normalize_source(p)
while numel(p)>1 && (p(end)=='/' || p(end)==char(92)),p(end)=[];end
if exist(fullfile(p,'rg_params.m'),'file')==2,return;end;error('W2C:Source','sourceRoot must directly contain rg_params.m.');
end
function local_restore(oldPath,oldPwd)
try,path(oldPath);catch,end;try,if exist(oldPwd,'dir')==7,cd(oldPwd);end;catch,end
end
function local_diary_off(),try,diary off;catch,end,end
function local_write_error(file,J)
fid=fopen(file,'w');if fid<0,return;end;c=onCleanup(@()fclose(fid));fprintf(fid,'Method: %s\nMode: %s\nIdentifier: %s\nMessage: %s\n\n%s',J.Method,J.Mode,J.ErrorId,J.ErrorMessage,J.ErrorStack); %#ok<NASGU>
end
function F=local_expand_failed(F,T)
fn=fieldnames(T);for k=1:numel(fn),if ~isfield(F,fn{k}),v=T.(fn{k});if ischar(v),F.(fn{k})='';elseif islogical(v),F.(fn{k})=false;else,F.(fn{k})=NaN;end;end;end
% Reorder to template field order for struct concatenation.
G=struct();for k=1:numel(fn),G.(fn{k})=F.(fn{k});end;F=G;
end
function local_write_readme(file,R)
fid=fopen(file,'w');if fid<0,return;end;c=onCleanup(@()fclose(fid));fprintf(fid,'Paper 1 -> Aerospace W2-C deterministic comparison\nGenerated: %s\nDecision: %s\nEvidence status: %s\n\n',R.Timestamp,R.Decision,R.EvidenceStatus);fprintf(fid,'Interpretation boundary: same certified rhoB=0.10 and same downstream loop; no Monte Carlo; no general controller ranking; development timing is not publication timing.\n'); %#ok<NASGU>
end
