function Report=run_p1_aerospace_w2d(sourceRoot,varargin)
%RUN_P1_AEROSPACE_W2D One-click predeclared targeted full-6DOF validation.
% MATLAB R2018a. Frozen Paper-1 source is read-only. No Monte Carlo is run.

p=inputParser;
p.addRequired('sourceRoot',@ischar);
p.addParameter('UseParallel',true,@(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('NumWorkers',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('SaveFullTrajectories',true,@(x)islogical(x)||ismember(x,[0 1]));
p.parse(sourceRoot,varargin{:});
A=p.Results;

packageRoot=fileparts(mfilename('fullpath'));
oldPath=path;
oldPwd=pwd;
cd(packageRoot);
addpath(sourceRoot,'-begin');
addpath(packageRoot,'-begin');
rehash path;
cleanupPath=onCleanup(@()local_restore(oldPath,oldPwd)); %#ok<NASGU>

sourceRoot=local_normalize_source(sourceRoot);
O=w2d_options();
O.UseParallel=logical(A.UseParallel);
O.RequestedWorkers=round(A.NumWorkers);
O.SaveFullTrajectories=logical(A.SaveFullTrajectories);

stamp=datestr(now,'yyyymmdd_HHMMSS');
outdir=[sourceRoot '_W2D_' stamp];
if exist(outdir,'dir')==7
    error('W2D:OutputExists','Output already exists: %s',outdir);
end
mkdir(outdir);
logfile=fullfile(outdir,'W2D_console_log.txt');
diary(logfile);
cleanupDiary=onCleanup(@()local_diary_off()); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('Paper 1 -> Aerospace | W2-D predeclared targeted constraint validation\n');
fprintf('============================================================\n');
fprintf('Version : %s\n',O.Version);
fprintf('MATLAB  : %s\n',version);
fprintf('Source  : %s\n',sourceRoot);
fprintf('Package : %s\n',packageRoot);
fprintf('Output  : %s\n',outdir);
fprintf('rhoB    : %.2f (in memory only)\n',O.RhoB);
fprintf('Pulse   : C2 quintic rise/fall, t0=%.2f s, ramp=%.2f s, hold=%.2f s\n',O.RampStart_s,O.RampDuration_s,O.PulseHold_s);
fprintf('Targets : T_min raw V=30 m/s; T_max raw V=46 m/s; N_plus raw gamma=+8 deg\n');
fprintf('No frozen source file will be modified. No Monte Carlo is run.\n');

Preflight=w2d_preflight(sourceRoot,outdir,O);
fprintf('PREFLIGHT PASS: exact frozen source, 18 hard rows, locked QP metric and target protocol.\n');

C0=cl6_params();
C0.R.rhoB=O.RhoB;
[WitnessRows,WitnessPass]=w2d_witness_reproduction(C0,O);
writetable(struct2table([WitnessRows{:}]),fullfile(outdir,'W2D_witness_reproduction.csv'));

fprintf('\nPART A -- W2-B trim-witness reproduction (algebraic, not closed-loop evidence)\n');
for k=1:numel(WitnessRows)
    r=WitnessRows{k};
    fprintf('  %-12s expected=%-7s scalarLimiter=%d qpSupport=%d QPpath=%s\n', ...
        r.Mode,r.ExpectedRow,r.ProposedExpectedLimiter,r.QPExpectedSupport,r.QPPath);
end
if ~WitnessPass
    error('W2D:Witness','Predeclared W2-B witness reproduction failed; do not run trajectories.');
end

methods=O.Methods;
modes=O.Modes;
N=numel(methods)*numel(modes);
jobMode=cell(N,1);
jobMethod=cell(N,1);
q=0;
for im=1:numel(modes)
    for it=1:numel(methods)
        q=q+1;
        jobMode{q}=modes{im};
        jobMethod{q}=methods{it};
    end
end

fprintf('\nPART B -- six paired targeted full-6DOF trajectories\n');
for j=1:N
    fprintf('  job %d: %-9s | %s\n',j,jobMethod{j},jobMode{j});
end

[P,pool]=w2d_parallel_prepare(O.UseParallel,O.RequestedWorkers,N,sourceRoot,packageRoot,outdir);
Jobs=cell(N,1);
try
    if P.Enabled
        parfor j=1:N
            Jobs{j}=w2d_run_job(jobMode{j},jobMethod{j},O);
        end
    else
        for j=1:N
            Jobs{j}=w2d_run_job(jobMode{j},jobMethod{j},O);
        end
    end
catch ME
    w2d_parallel_finish(pool,P,sourceRoot,packageRoot,outdir);
    pool=[];
    rethrow(ME);
end
PoolFinish=w2d_parallel_finish(pool,P,sourceRoot,packageRoot,outdir);
pool=[];
fprintf('%s\n',PoolFinish.Message);

Cases=struct();
Entries=cell(N,1);
allOK=true;
for j=1:N
    J=Jobs{j};
    Entries{j}=J.Entry;
    key=[J.Method '_' J.Mode];
    key=strrep(key,'-','_');
    if J.Success
        Cases.(key)=J.S;
        E=J.Entry;
        fprintf('[%-9s | %-12s] COMPLETED: RMS=%.5f m peak=%.5f m intervention=%.2f%% expected=%s support=%d limiter=%d maxVB/rho=%.9f\n', ...
            J.Method,J.Mode,E.PosRMS_Full_m,E.PosPeak_Full_m,E.GovInterventionPct, ...
            E.ExpectedRow,E.ExpectedQPSupportSamples,E.ExpectedScalarLimiterSamples,E.MaxVB_over_rhoB);
        if O.SaveFullTrajectories
            S=J.S; %#ok<NASGU>
            save(fullfile(outdir,[key '.mat']),'S','-v7.3');
        end
    else
        allOK=false;
        fprintf(2,'[%-9s | %-12s] FAILED: %s -- %s\n',J.Method,J.Mode,J.ErrorId,J.ErrorMessage);
        local_write_error(fullfile(outdir,[key '_ERROR.txt']),J);
    end
end

Template=[];
for j=1:N
    if Jobs{j}.Success
        Template=Jobs{j}.Entry;
        break;
    end
end
for j=1:N
    if ~Jobs{j}.Success && ~isempty(Template)
        Entries{j}=local_expand_failed(Jobs{j}.Entry,Template);
    end
end
writetable(struct2table([Entries{:}]),fullfile(outdir,'W2D_case_summary.csv'));

PairRows=cell(0,1);
RoleRows=cell(0,1);
PhiRows=cell(0,1);
MarginRows=cell(0,1);
InitialRows=cell(0,1);

if allOK
    fprintf('\nPART C -- own-trajectory targeted evidence audit\n');
    for im=1:numel(modes)
        SP=Cases.(['Proposed_' modes{im}]);
        SQ=Cases.(['MatchedQP_' modes{im}]);
        EP=w2d_summarize_case(SP,O);
        EQ=w2d_summarize_case(SQ,O);
        pr=w2d_pairwise_row(SP,SQ,EP,EQ);
        PairRows{end+1}=pr; %#ok<AGROW>
        fprintf('%-12s expected=%-7s Proposed limiter=%d | QP support=%d activeSet=%d | dRMS=%+.5f m dPeak=%+.5f m\n', ...
            modes{im},pr.ExpectedRow,pr.ProposedExpectedLimiterSamples,pr.QPExpectedSupportSamples, ...
            pr.QPExpectedActiveSetSamples,pr.QP_minus_Proposed_RMS_m,pr.QP_minus_Proposed_Peak_m);
    end

    for j=1:N
        S=Cases.([Jobs{j}.Method '_' Jobs{j}.Mode]);
        rr=w2d_role_rows(S,O);
        RoleRows=[RoleRows;rr(:)]; %#ok<AGROW>
        [pp,mm]=w2d_margin_rows(S);
        PhiRows=[PhiRows;pp(:)]; %#ok<AGROW>
        MarginRows=[MarginRows;mm(:)]; %#ok<AGROW>
        InitialRows{end+1}=struct( ... %#ok<AGROW>
            'Method',S.Method,'Mode',S.mode,'VB0_over_rhoB',S.I0.VB0_over_rho, ...
            'CertificateCompatible',S.I0.VB0<=S.C.R.rhoB+O.InitialBackupTol);
    end

    writetable(struct2table([PairRows{:}]),fullfile(outdir,'W2D_pairwise_comparison.csv'));
    writetable(struct2table([RoleRows{:}]),fullfile(outdir,'W2D_constraint_roles.csv'));
    writetable(struct2table([PhiRows{:}]),fullfile(outdir,'W2D_affine_row_minima.csv'));
    writetable(struct2table([MarginRows{:}]),fullfile(outdir,'W2D_reference_margin_minima.csv'));
    writetable(struct2table([InitialRows{:}]),fullfile(outdir,'W2D_initial_compatibility.csv'));
else
    fprintf(2,'PART C skipped because at least one targeted trajectory failed.\n');
end

CorePass=allOK && WitnessPass;
if allOK
    for j=1:N
        E=Jobs{j}.Entry;
        CorePass=CorePass && E.ReferenceGatePass && E.ActuatorBoxPass && ...
            E.QPKKTGatePass && E.InitialVB_over_rhoB<=1+O.InitialBackupTol;
    end

    support=zeros(numel(PairRows),1);
    for k=1:numel(PairRows)
        support(k)=PairRows{k}.QPExpectedSupportSamples;
    end
    if all(support>0)
        EvidenceStatus='FULL_PREDECLARED_TARGET_SUPPORT_CONFIRMED';
    elseif any(support>0)
        EvidenceStatus='PARTIAL_PREDECLARED_TARGET_SUPPORT_CONFIRMED';
    else
        EvidenceStatus='NO_PREDECLARED_TARGET_SUPPORT_FOUND';
    end
else
    EvidenceStatus='NOT_EVALUATED_DUE_TO_TRAJECTORY_FAILURE';
end

if CorePass
    Decision='PASS_W2D_TARGETED_CORE';
else
    Decision='HOLD_W2D_CORE_FAILURE';
end

if CorePass && ~strcmp(EvidenceStatus,'NO_PREDECLARED_TARGET_SUPPORT_FOUND')
    NextStatus='READY_TO_FREEZE_EVIDENCE_AND_DESIGN_PAIRED_MONTE_CARLO';
elseif CorePass
    NextStatus='CLAIM_SHRINK_DECISION_REQUIRED_BEFORE_MONTE_CARLO';
else
    NextStatus='DO_NOT_ENTER_MONTE_CARLO';
end

fprintf('\n============================================================\n');
fprintf('W2-D DECISION: %s\n',Decision);
fprintf('EVIDENCE STATUS: %s\n',EvidenceStatus);
fprintf('NEXT STATUS: %s\n',NextStatus);
fprintf('============================================================\n');
fprintf('Target cases were predeclared from W2-B finite witnesses; no post-result retuning is performed by this runner.\n');

Report=struct('Timestamp',datestr(now,31),'Version',O.Version,'SourceRoot',sourceRoot, ...
    'PackageRoot',packageRoot,'OutputDir',outdir,'Options',O,'Preflight',Preflight, ...
    'WitnessRows',{WitnessRows},'WitnessPass',WitnessPass,'Parallel',P,'ParallelFinish',PoolFinish, ...
    'Jobs',{Jobs},'Cases',Cases,'PairRows',{PairRows},'RoleRows',{RoleRows}, ...
    'PhiRows',{PhiRows},'MarginRows',{MarginRows},'Decision',Decision, ...
    'EvidenceStatus',EvidenceStatus,'NextStatus',NextStatus);
save(fullfile(outdir,'W2D_Report.mat'),'Report','-v7.3');
local_write_readme(fullfile(outdir,'W2D_README_RESULTS.txt'),Report);
diary off;

zipPath=[outdir '.zip'];
try
    zip(zipPath,{'*'},outdir);
    fprintf('Output ZIP: %s\n',zipPath);
catch ME
    fprintf(2,'ZIP creation failed: %s\n',ME.message);
end
fprintf('Output directory: %s\n',outdir);
end

function p=local_normalize_source(p)
while numel(p)>1 && (p(end)=='/' || p(end)==char(92))
    p(end)=[];
end
if exist(fullfile(p,'rg_params.m'),'file')==2
    return;
end
error('W2D:Source','sourceRoot must directly contain rg_params.m.');
end

function local_restore(oldPath,oldPwd)
try
    path(oldPath);
catch
end
try
    if exist(oldPwd,'dir')==7
        cd(oldPwd);
    end
catch
end
end

function local_diary_off()
try
    diary off;
catch
end
end

function local_write_error(file,J)
fid=fopen(file,'w');
if fid<0
    return;
end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Method: %s\nMode: %s\nIdentifier: %s\nMessage: %s\n\n%s', ...
    J.Method,J.Mode,J.ErrorId,J.ErrorMessage,J.ErrorStack);
end

function F=local_expand_failed(F,T)
fn=fieldnames(T);
for k=1:numel(fn)
    if ~isfield(F,fn{k})
        v=T.(fn{k});
        if ischar(v)
            F.(fn{k})='';
        elseif islogical(v)
            F.(fn{k})=false;
        else
            F.(fn{k})=NaN;
        end
    end
end
G=struct();
for k=1:numel(fn)
    G.(fn{k})=F.(fn{k});
end
F=G;
end

function local_write_readme(file,R)
fid=fopen(file,'w');
if fid<0
    return;
end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Paper 1 -> Aerospace W2-D targeted deterministic validation\n');
fprintf(fid,'Generated: %s\nDecision: %s\nEvidence: %s\nNext: %s\n\n', ...
    R.Timestamp,R.Decision,R.EvidenceStatus,R.NextStatus);
fprintf(fid,'The targets are predeclared from W2-B finite witnesses. This is closed-loop deterministic evidence, not Monte Carlo evidence or a general controller ranking.\n');
end
