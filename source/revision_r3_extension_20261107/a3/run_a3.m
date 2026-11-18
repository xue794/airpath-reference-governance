function run_a3()
h=fileparts(mfilename('fullpath')); root=fileparts(fileparts(h));
core=fullfile(root,'FWUAV','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel');
pkg=fullfile(root,'FWUAV','P1_Aerospace_W2E_R2018a_v1_0'); aoa=fullfile(root,'contribution_probe_20261106','aoa');
addpath(core,pkg,aoa,h); assert(strcmp(version('-release'),'2018a'));
O=w2e_options(); P=jsondecode(fileread(fullfile(h,'protocol.json')));
ids=P.ids(:).'; arms={'N10','D16','Am'}; steps=[.01,.005];
assert(isequal(ids,[22 32 33]));
E=readtable(fullfile(fileparts(h),'a2','event_comparison.csv'));
jobs=struct('id',{},'arm',{},'step',{});
for id=ids
    tag=sprintf('id_%03d.mat',id); b=load(fullfile(fileparts(h),'a1','runs',tag),'S');
    a=load(fullfile(fileparts(h),'a2','runs',tag),'u16','u18'); S=b.S; C=S.C;
    [x0,~]=cl6_initial_state(O.Mode,C);
    assert(isequal(x0,S.x(1,:).') && isequal(S.t,(0:.05:120).') && C.MaxStep==.01 && ...
       C.RelTol==1e-8 && C.AbsTol==1e-10 && C.R.rhoB==.10 && C.R.DeltaN==10 && ...
       C.R.Tmin+C.R.DeltaT==3 && C.R.Tmax-C.R.DeltaT==57 && strcmp(S.mode,O.Mode));
    assert(isequal(w2e_rhs(0,x0,O.Mode,C,'MatchedQP',O),a3_rhs(0,x0,O.Mode,C,'N10',O)));
    ix=unique([1;E{E{:,1}==id,2}]);
    for j=1:numel(ix)
        k=ix(j); x=S.x(k,:).'; n=a3_snapshot(S.t(k),x,O.Mode,C,'N10',O);
        d=a3_snapshot(S.t(k),x,O.Mode,C,'D16',O); m=a3_snapshot(S.t(k),x,O.Mode,C,'Am',O);
        assert(isequal(n.Gov.u,a.u18(k,:).') && isequal(d.Gov.u,a.u16(k,:).'));
        common=[1:13,16:18];assert(isequal(d.Gov.appliedG,m.Gov.appliedG(common,:)) && ...
          isequal(d.Gov.appliedC,m.Gov.appliedC(common)) && numel(m.Gov.phi)==18);
    end
    for si=1:2
        for ai=1:3
            jobs(end+1)=struct('id',id,'arm',arms{ai},'step',steps(si)); %#ok<AGROW>
        end
    end
end
fprintf('A3 preflight passed archived sample, RHS, force reserve, row ordering and A2 checks.\n');
[info,pool]=w2e_parallel_prepare(true,numel(jobs),numel(jobs),core,pkg,fullfile(h,'identity'));
assert(info.Enabled,'A3:Pool','Parallel pool unavailable');
parfor j=1:numel(jobs)
    addpath(h,aoa); a3_run_one(h,jobs(j),O);
end
w2e_parallel_finish(pool,info,core,pkg,fullfile(h,'identity'));
fprintf('A3 batch returned. Inspect all job status files.\n');
end
