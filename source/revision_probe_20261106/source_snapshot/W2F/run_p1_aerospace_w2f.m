function Report = run_p1_aerospace_w2f(sourceRoot)
%RUN_P1_AEROSPACE_W2F Publication-grade serial timing benchmark, MATLAB R2018a.
%
% Report=run_p1_aerospace_w2f(...frozen Paper-1 source root...)
%
% This benchmark intentionally DOES NOT use the 56-worker pool. Runtime claims
% must be measured on the MATLAB client in isolation. No ODE or Monte Carlo runs.
packageRoot=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(sourceRoot),sourceRoot=local_find_source();end
if isa(sourceRoot,'string') && isscalar(sourceRoot),sourceRoot=char(sourceRoot);end
assert(ischar(sourceRoot) && size(sourceRoot,1)==1,'W2F:Path','sourceRoot must be a character path.');
while numel(sourceRoot)>1 && (sourceRoot(end)=='/' || sourceRoot(end)==char(92)),sourceRoot=sourceRoot(1:end-1);end
if exist(fullfile(sourceRoot,'rg_params.m'),'file')~=2
    nested=fullfile(sourceRoot,'T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel');
    if exist(fullfile(nested,'rg_params.m'),'file')==2,sourceRoot=nested;else,error('W2F:SourceRoot','Select folder containing rg_params.m.');end
end
oldPath=path;pathCleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(sourceRoot,'-begin');addpath(packageRoot,'-begin');rehash path;
clear rg_params rg_constraint_affine rg_scalar_governor rg_backup_control rg_backup_state rg_nominal_control;
clear rg_angle_error rg_N_map rg_nominal_aero rg_alpha_root t2_gna_params t2_isa_density rg_solve_lyap;
[parent,name]=fileparts(sourceRoot);if isempty(parent),parent=pwd;end
stamp=datestr(now,'yyyymmdd_HHMMSS');outdir=fullfile(parent,[name '_W2F_PublicationTiming_rho010_v1_0_' stamp]);
[ok,msg]=mkdir(outdir);assert(ok,'W2F:Output','Cannot create output: %s',msg);
diary(fullfile(outdir,'W2F_console_log.txt'));diaryCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
O=w2f_options();
Report=struct('Version',O.Version,'Protocol',O.Protocol,'Timestamp',datestr(now,31),'Environment',version, ...
    'SourceRoot',sourceRoot,'PackageRoot',packageRoot,'OutputDir',outdir,'OutputZip',[outdir '.zip'], ...
    'Options',O,'Machine',w2f_machine_info(),'ThreadControl',[],'Preflight',[],'Benchmark',[], ...
    'Decision','RUNNING','ErrorId','','ErrorMessage','');
oldThreads=[];threadApplied=false;
fprintf('\n============================================================\n');
fprintf('Paper 1 -> Aerospace | W2-F publication computational timing\n');
fprintf('============================================================\n');
fprintf('Version  : %s\nProtocol : %s\nMATLAB   : %s\nSource   : %s\nOutput   : %s\n',O.Version,O.Protocol,version,sourceRoot,outdir);
fprintf('Mode     : CLIENT-SERIAL ONLY; no ODE, no Monte Carlo, no parfor.\n');
fprintf('Bank     : frozen W2-B 425-state common-state benchmark.\n');
fprintf('rhoB     : 0.10 in memory only.\n');
try
    if O.RequireNoParallelPool && exist('gcp','file')==2
        p=gcp('nocreate');
        if ~isempty(p)
            error('W2F:ParallelPoolOpen','A parallel pool is open (%d workers). Close it before publication timing.',p.NumWorkers);
        end
    end
    if O.ForceSingleComputeThread
        try
            oldThreads=maxNumCompThreads;
            maxNumCompThreads(1);
            threadApplied=true;
            fprintf('Thread gate: MATLAB client computation capped at one thread for timed region (previous max=%g).\n',oldThreads);
        catch MEthread
            error('W2F:ThreadControl','Could not enforce one MATLAB computation thread: %s',MEthread.message);
        end
    end
    threadCleanup=onCleanup(@()local_restore_threads(threadApplied,oldThreads)); %#ok<NASGU>
    fprintf('\nPART A -- exact source/bank/interface/correctness gate\n');
    [R,Bank,Gate]=w2f_preflight(sourceRoot,outdir,O);Report.Preflight=Gate;
    fprintf('PASS: 425 frozen states; selector cores reproduce source implementations; 18 rows; rhoB=0.10.\n');
    fprintf('QP paths and case composition written before timing.\n');
    save(fullfile(outdir,'W2F_Report.mat'),'Report','-v7');
    fprintf('\nPART B -- warmed, blocked, serial publication timing\n');
    B=w2f_benchmark(R,Bank,O,outdir);Report.Benchmark=B;
    Report.Decision='PASS_W2F_PUBLICATION_TIMING_COMPLETED';
catch ME
    Report.Decision='HOLD_W2F_TIMING_ERROR';Report.ErrorId=ME.identifier;Report.ErrorMessage=ME.message;
    fprintf(2,'\n%s: %s\n',ME.identifier,ME.message);
    for i=1:numel(ME.stack),fprintf(2,'  %s line %d\n',ME.stack(i).name,ME.stack(i).line);end
end
local_restore_threads(threadApplied,oldThreads);threadApplied=false;
save(fullfile(outdir,'W2F_Report.mat'),'Report','-v7');
local_environment_text(Report,outdir);
local_summary(Report,outdir);
fprintf('\n============================================================\nW2-F DECISION: %s\n',Report.Decision);
if strcmp(Report.Decision,'PASS_W2F_PUBLICATION_TIMING_COMPLETED')
    P=Report.Benchmark.PublicationTable;
    fprintf('Publication medians: selector Scalar=%.6f ms, QP=%.6f ms, ratio=%.3f\n',P.ScalarSelectorMedian_ms,P.QPSelectorMedian_ms,P.QP_over_ScalarSelectorMedian);
    fprintf('                     whole    Scalar=%.6f ms, QP=%.6f ms, ratio=%.3f\n',P.ScalarWholeMedian_ms,P.QPWholeMedian_ms,P.QP_over_ScalarWholeMedian);
end
fprintf('Output directory: %s\n',outdir);
diary off;
try
    L=dir(outdir);names={L(~[L.isdir]).name};zip(Report.OutputZip,names,outdir);fprintf('Output ZIP: %s\n',Report.OutputZip);
catch MEzip
    warning('W2F:Zip','Could not ZIP outputs: %s',MEzip.message);
end
end

function local_restore_threads(applied,oldThreads)
if applied && ~isempty(oldThreads)
    try,maxNumCompThreads(oldThreads);catch,end
end
end
function p=local_find_source()
known='C:\FWUAV\T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel';
if exist(fullfile(known,'rg_params.m'),'file')==2,p=known;return;end
q=which('rg_params');if ~isempty(q),p=fileparts(q);return;end
p=uigetdir(pwd,'Select frozen Paper 1 source folder containing rg_params.m');
if isequal(p,0),error('W2F:Cancelled','Source selection cancelled.');end
end
function local_environment_text(R,outdir)
fid=fopen(fullfile(outdir,'W2F_environment.txt'),'w');if fid<0,return;end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Protocol: %s\nVersion: %s\nTimestamp: %s\nMATLAB: %s\nComputer: %s\n',R.Protocol,R.Version,R.Timestamp,R.Machine.MATLABVersion,R.Machine.Computer);
fprintf(fid,'Hostname: %s\nProcessorIdentifier: %s\nNUMBER_OF_PROCESSORS: %s\nfeature(numcores): %.0f\n',R.Machine.Hostname,R.Machine.ProcessorIdentifier,R.Machine.NumberOfProcessors,R.Machine.NumCoresFeature);
fprintf(fid,'OS: %s\nCPU query:\n%s\n',R.Machine.OS,R.Machine.CPUQuery);
fprintf(fid,'Timed region requested one MATLAB compute thread and no parallel pool.\n');
end
function local_summary(R,outdir)
fid=fopen(fullfile(outdir,'W2F_SUMMARY.txt'),'w');if fid<0,return;end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Paper 1 W2-F Publication Timing\nDecision: %s\nProtocol: %s\n',R.Decision,R.Protocol);
fprintf(fid,'Frozen common-state bank: 425 W2-B snapshots. No new 6DOF/MC runs.\n');
fprintf(fid,'Client serial, no parallel pool, one MATLAB compute thread for timed region.\n');
fprintf(fid,'Warm-up=%d passes; measured=%d blocks x %d passes/block.\n',R.Options.WarmupPasses,R.Options.NumBlocks,R.Options.PassesPerBlock);
fprintf(fid,'Each method/timed item therefore has %d measured calls.\n',425*R.Options.NumBlocks*R.Options.PassesPerBlock);
if strcmp(R.Decision,'PASS_W2F_PUBLICATION_TIMING_COMPLETED')
    P=R.Benchmark.PublicationTable;
    fprintf(fid,'Constraint build median: %.9f ms/call\n',P.ConstraintBuildMedian_ms);
    fprintf(fid,'Backup control median: %.9f ms/call\n',P.BackupControlMedian_ms);
    fprintf(fid,'Scalar selector median: %.9f ms/call\n',P.ScalarSelectorMedian_ms);
    fprintf(fid,'QP selector median: %.9f ms/call\n',P.QPSelectorMedian_ms);
    fprintf(fid,'QP/scalar selector ratio median: %.6f\n',P.QP_over_ScalarSelectorMedian);
    fprintf(fid,'Scalar whole governor median: %.9f ms/call\n',P.ScalarWholeMedian_ms);
    fprintf(fid,'QP whole governor median: %.9f ms/call\n',P.QPWholeMedian_ms);
    fprintf(fid,'QP/scalar whole-governor ratio median: %.6f\n',P.QP_over_ScalarWholeMedian);
end
fprintf(fid,'Report medians and IQR/block dispersion; do not use parallel MC wall time as algorithm timing.\n');
if ~isempty(R.ErrorId),fprintf(fid,'ERROR: %s -- %s\n',R.ErrorId,R.ErrorMessage);end
end
