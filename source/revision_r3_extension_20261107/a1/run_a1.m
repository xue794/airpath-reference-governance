function run_a1()
here=fileparts(mfilename('fullpath')); project=fileparts(fileparts(here));
core=fullfile(project,'FWUAV','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel');
pkg=fullfile(project,'FWUAV','P1_Aerospace_W2E_R2018a_v1_0');
archive=[core '_W2E_PairedMC500_rho010_v1_0'];
for dir={'identity','runs','checks','logs'}
    folder=fullfile(here,dir{1}); if ~isdir(folder),mkdir(folder);end
end
addpath(core,'-begin');addpath(pkg,'-begin');rehash path;
assert(strcmp(version('-release'),'2018a'),'A1:Version','Requires MATLAB R2018a');
O=w2e_options(); P=w2e_preflight(core,fullfile(here,'identity'),O);
assert(strcmpi(P.SampleHash,'858b3a0fcaa70f5b76791cc400991d5ff8e02114f0aff1005a2568461fc69b46'));
ids=[22 32 33 80 86 91 108 113 164 184 236 322 374 386 397 472 483];
F=load(fullfile(archive,'W2E_CHECKPOINT.mat'),'PairResultsSaved','ProtocolIDSaved','SampleHashSaved');
assert(strcmpi(F.SampleHashSaved,P.SampleHash) && strcmp(F.ProtocolIDSaved,O.ProtocolID));
Eold=cellfun(@(r)r.MatchedQP,F.PairResultsSaved,'UniformOutput',false);
chosen=find(cellfun(@(e)e.CorePass && e.RoleCount(15)>0,Eold));
assert(isequal(chosen(:),ids(:)) && sum(cellfun(@(e)e.RoleCount(15),Eold(ids)))==21);
assert(all(cellfun(@(e)isequal(e.BindingCount(15),e.RoleCount(15)),Eold(ids))));
Samples=P.Samples(ids); expected=Eold(ids); Cbase=P.Cbase; Rmc=P.Rmc;
save(fullfile(here,'identity','frozen.mat'),'Samples','expected','Cbase','Rmc','O','ids','-v7');
fprintf('PREFLIGHT PASS %d source files, %d frozen samples, 17 selected, 21 supports.\n',numel(P.SourceRows),numel(P.Samples));
pending=find(~arrayfun(@(id)exist(fullfile(here,'checks',sprintf('id_%03d.json',id)),'file')==2,ids));
if isempty(pending),fprintf('All cases have recorded outcomes.\n');return;end
[info,pool]=w2e_parallel_prepare(true,numel(pending),numel(pending),core,pkg,fullfile(here,'identity'));
if info.Enabled
    parfor j=1:numel(pending)
        k=pending(j);run_one(here,Samples(k),expected{k},Cbase,Rmc,O);
    end
else
    for k=pending(:).',run_one(here,Samples(k),expected{k},Cbase,Rmc,O);end
end
w2e_parallel_finish(pool,info,core,pkg,fullfile(here,'identity'));
fprintf('A1 RUNS FINISHED\n');
end

function run_one(here,sample,old,Cbase,Rmc,O)
id=sample.index; tag=sprintf('id_%03d',id);
logfile=fullfile(here,'logs',[tag '.txt']);
fid=fopen(logfile,'a'); cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s START ID %d\n',datestr(now,31),id);
try
    C=rev41_apply_sample(Cbase,sample,Rmc);
    assert(isequaln(C.R,Cbase.R) && isequaln(C.A,Cbase.A) && C.R.rhoB==.10);
    S=w2e_simulate_trial('MatchedQP',C,O); E=w2e_summarize_run(S,sample,Rmc,O);
    % Persist the original output before the separate algebraic export.
    raw=fullfile(here,'runs',[tag '.mat']); partial=[raw '.partial.mat'];
    save(partial,'S','sample','E','old','-v7');movefile(partial,raw,'f');
    fprintf(fid,'%s ORIGINAL OUTPUT SAVED\n',datestr(now,31));
    n=numel(S.t); [x0,~]=cl6_initial_state(O.Mode,C);
    samegrid=isequal(S.t,(0:.05:120).') && isequal(size(S.x),[2401 36]);
    sameinit=isequal(S.x(1,:).',x0); names=w2b_row_names();
    T=struct('Index',id,'N',n,'Status',E.Status,'StatusMatch',strcmp(E.Status,old.Status), ...
        'CorePass',E.CorePass,'CoreMatch',isequal(E.CorePass,old.CorePass), ...
        'CompletedMatch',isequal(E.Completed,old.Completed),'GridMatch',samegrid,'InitialMatch',sameinit, ...
        'NamesMatch',isequal(S.Y.constraintNames(:),names(:)), ...
        'RoleMatch',isequal(E.RoleCount,old.RoleCount),'BindingMatch',isequal(E.BindingCount,old.BindingCount), ...
        'Support',E.RoleCount(15),'Binding',E.BindingCount(15), ...
        'RMSFull',E.PosRMSFull_m,'OldFull',old.PosRMSFull_m,'RMSLate',E.PosRMSLate_m,'OldLate',old.PosRMSLate_m, ...
        'KKT',max([E.QPMaxPrimal,E.QPMaxStationarity,E.QPMaxDualViolation,E.QPMaxComplementarity]), ...
        'RowViolation',E.MaxNormalizedPhiViolation,'MaxVBoverRho',E.MaxVB_over_rhoB, ...
        'Integration_s',S.IntegrationWall_s,'Output_s',S.OutputWall_s);
    T.RMSPass=abs(T.RMSFull-T.OldFull)<=max(1e-6,1e-8*abs(T.OldFull)) && ...
        abs(T.RMSLate-T.OldLate)<=max(1e-6,1e-8*abs(T.OldLate));
    assert(isequal(names(14:15),{'N_minus';'N_plus'}) && strcmp(names{18},'backup'));
    affineC=zeros(n,18); affineG=zeros(n,18,3); dual=zeros(n,18); dualScaled=zeros(n,18);
    VBoverRho=zeros(n,1); alphaRef=zeros(n,1); alphaFeedback=zeros(n,1); gapRef=zeros(n,2); gapFeedback=zeros(n,2);
    maxPhi=0; maxU=0; maxReq=0; flags=true;
    for k=1:n
        snap=w2e_control_snapshot(S.t(k),S.x(k,:).',O.Mode,C,'MatchedQP',O);
        ac=snap.Gov.C; affineC(k,:)=ac.c.';affineG(k,:,:)=reshape(ac.G,[1,18,3]);
        dual(k,:)=snap.Gov.Solver.Dual.';dualScaled(k,:)=snap.Gov.Solver.DualScaled.';
        supportMask=false(18,1);supportMask(snap.Gov.Solver.SupportRows)=true;
        flags=flags && isequal(ac.names(:),names(:)) && isequal(S.Y.qpSupport(k,:).',supportMask);
        maxPhi=max(maxPhi,norm(ac.c+ac.G*S.Y.uGov(k,:).'-S.Y.phi(k,:).',inf)/(1+norm(S.Y.phi(k,:),inf)));
        maxU=max(maxU,norm((snap.Gov.u-S.Y.uGov(k,:).')./O.QP.Scale,inf));
        maxReq=max(maxReq,norm((snap.u0-S.Y.u0(k,:).')./O.QP.Scale,inf));
        xi=rg_backup_state(S.Y.xf(k,:).',C.R);VBoverRho(k)=xi.'*C.R.PB*xi/C.R.rhoB;
        alphaRef(k)=ac.alpha; alphaFeedback(k)=snap.Cmd.alpha;
        gapRef(k,:)=[ac.hNminus,ac.hNplus]+C.R.DeltaN;
        lo=rg_N_map(C.R.alphaLo,snap.Air.V,snap.Cmd.X,C.R);hi=rg_N_map(C.R.alphaHi,snap.Air.V,snap.Cmd.X,C.R);
        gapFeedback(k,:)=[snap.Cmd.N-lo.value,hi.value-snap.Cmd.N];
    end
    T.SnapshotConsistent=flags && maxPhi<=1e-7 && maxU<=1e-8 && maxReq<=1e-8;
    T.MaxPhiReplay=maxPhi;T.MaxUReplay=maxU;T.MaxRequestReplay=maxReq;
    T.Pass=all([T.StatusMatch,T.CoreMatch,T.CorePass,T.CompletedMatch,E.Completed,samegrid,sameinit,T.NamesMatch, ...
        T.RoleMatch,T.BindingMatch,T.RMSPass,T.KKT<=O.KKTGateTol,T.RowViolation<=O.NormalizedPhiTol, ...
        T.MaxVBoverRho<=1+O.BackupTol,T.SnapshotConsistent]);
    partial=fullfile(here,'runs',[tag '_rows.partial.mat']); out=fullfile(here,'runs',[tag '_rows.mat']);
    save(partial,'affineC','affineG','names','dual','dualScaled','VBoverRho','alphaRef','alphaFeedback','gapRef','gapFeedback','T','-v7');
    movefile(partial,out,'f');
catch ME
    T=struct('Index',id,'Pass',false,'ErrorId',ME.identifier,'Error',getReport(ME,'extended','hyperlinks','off'));
end
out=fullfile(here,'checks',[tag '.json']); f=fopen([out '.tmp'],'w');fwrite(f,jsonencode(T),'char');fclose(f);movefile([out '.tmp'],out,'f');
fprintf(fid,'%s END pass=%d\n',datestr(now,31),T.Pass);
end
