function run_a2()
% A2 holds the original 36-state trajectory fixed; it performs no integration.
here=fileparts(mfilename('fullpath')); stage=fileparts(here); root=fileparts(stage);
a1=fullfile(stage,'a1'); core=fullfile(root,'FWUAV','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel');
pkg=fullfile(root,'FWUAV','P1_Aerospace_W2E_R2018a_v1_0'); b=fullfile(root,'FWUAV','P1_Aerospace_W2B_R2018a_v1_1');
assert(strcmp(version('-release'),'2018a')); addpath(b,'-end');addpath(core,'-begin');addpath(pkg,'-begin');
assert(strcmpi(which('w2b_admission'),fullfile(b,'w2b_admission.m')));
if ~isdir(fullfile(here,'runs')),mkdir(fullfile(here,'runs'));end
if ~isdir(fullfile(here,'identity')),mkdir(fullfile(here,'identity'));end
L=load(fullfile(a1,'identity','frozen.mat'),'ids','O');
live=w2e_options();assert(isequaln(rmfield(L.O,'SampleSource'),rmfield(live,'SampleSource')));
assert(strcmpi(w2e_sha256(live.SampleSource),'858b3a0fcaa70f5b76791cc400991d5ff8e02114f0aff1005a2568461fc69b46'));
ids=L.ids;
[info,pool]=w2e_parallel_prepare(true,8,numel(ids),core,pkg,fullfile(here,'identity'));
if info.Enabled
    parfor j=1:numel(ids), replay_case(a1,here,ids(j),L.O);end
else
    for j=1:numel(ids), replay_case(a1,here,ids(j),L.O);end
end
w2e_parallel_finish(pool,info,core,pkg,fullfile(here,'identity'));
fprintf('A2_REPLAY_COMPLETE\n');
end

function replay_case(a1,here,id,O)
tag=sprintf('id_%03d',id);
L=load(fullfile(a1,'runs',[tag '.mat']),'S'); S=L.S;
R=load(fullfile(a1,'runs',[tag '_rows.mat']));
keep=[1:13 16:18]; ends=[14 15]; Q=O.QP; scale=Q.Scale; n=numel(S.t);
assert(R.T.Pass && S.C.MC.sampleIndex==id && S.C.R.rhoB==.1);
assert(n==2401 && isequal(S.t,(0:.05:120).') && isequal(R.names,w2b_row_names()));
assert(isequal(R.names(14:18),{'N_minus';'N_plus';'T_min';'T_max';'backup'}));
% One column-major numeric vector per sample, with separate three-input arrays.
fields={'u0_match','u18_match','phi_match','coeff_match','support_match','chart_ok','admission_ok',...
    'full_norm_viol','kept_norm_viol','kkt18','kkt16','objective18','objective16',...
    'd','phi16_lower','phi16_upper','norm_lower','norm_upper','tol_lower','tol_upper',...
    'support_upper','support_backup','retained_backup_support','replay_ok','base_ok','candidate',...
    'check_enum','enum_kkt18','enum_kkt16','enum_d18','enum_d16','enum_verified','accepted'};
M=nan(n,numel(fields));u0=nan(n,3);u18=nan(n,3);u16=nan(n,3);
deletion_dual=nan(n,16);enum_u18=nan(n,3);enum_u16=nan(n,3);status=repmat({'not_evaluated'},n,1);
for k=1:n
    try
        sn=w2e_control_snapshot(S.t(k),S.x(k,:).',O.Mode,S.C,'MatchedQP',O);
        F=sn.Gov; A=F.C; original=F.Solver; u=F.u; un=sn.u0;
        H=w2b_admission(S.Y.xf(k,:).',S.C.R,A,Q);
        support=false(18,1);support(original.SupportRows)=true;
        du0=norm((un-S.Y.u0(k,:).')./scale,inf);
        du18=norm((u-S.Y.uGov(k,:).')./scale,inf);
        phi_match=norm(F.phi-S.Y.phi(k,:).',inf)/(1+norm(S.Y.phi(k,:),inf));
        oldG=reshape(R.affineG(k,:,:),18,3);
        coeff=isequal(A.c,R.affineC(k,:).') && isequal(A.G,oldG) && isequal(A.names,R.names);
        sup=isequal(support,S.Y.qpSupport(k,:).');
        fullviol=viol(A,u,1:18);k18=kkt(original);
        replay=original.Success && du0<=Q.ReferenceMatchTol && du18<=Q.ReferenceMatchTol && ...
            phi_match<=1e-7 && coeff && sup && sn.chart.regular && H.Compatible && ...
            H.VB_over_rhoB<=1+O.BackupTol && fullviol<=O.NormalizedPhiTol && k18<=O.KKTGateTol;
        M(k,1:8)=[du0,du18,phi_match,coeff,sup,sn.chart.regular,H.Compatible,fullviol];
        M(k,10)=k18;M(k,21:22)=support([15 18]).';M(k,24)=replay;
        u0(k,:)=un.';u18(k,:)=u.';
        if ~replay,status{k}='full_replay_gate';continue;end
        opt=Q;opt.FeasibleStart=F.uB;
        [v,D]=w2b_qp3(un,-A.G(keep,:),A.c(keep),scale,opt);
        if ~D.Success,status{k}=['subset_' D.Status];continue;end
        u16(k,:)=v.';deletion_dual(k,:)=D.Dual.';
        change=norm((v-u)./scale,inf); raw=A.c(ends)+A.G(ends,:)*v;
        normneg=max(0,-raw./(1+abs(A.c(ends))+abs(A.G(ends,:))*abs(v)));
        rt=raw_tolerance(A,ends,un,v,Q);
        kept=viol(A,v,keep);k16=kkt(D);
        nested=D.Objective<=original.Objective+Q.ObjectiveTol*(1+original.Objective);
        base=replay && kept<=O.NormalizedPhiTol && k16<=O.KKTGateTol && nested;
        changed=change>Q.ReferenceMatchTol;
        normalized=normneg>O.NormalizedPhiTol; outside=-raw>rt;
        candidate=base && changed && any(normalized & outside);
        ratios=[change/Q.ReferenceMatchTol;normneg/O.NormalizedPhiTol;-raw./max(rt,realmin)];
        near=any(ratios>=.5 & ratios<=2);
        check=any(support(ends)) || changed || any(normalized | outside) || near || ~base;
        M(k,9)=kept;M(k,11:20)=[k16,original.Objective,D.Objective,change,raw.',normneg.',rt.'];
        M(k,23)=ismember(16,D.SupportRows);M(k,25:27)=[base,candidate,check];M(k,33)=0;
        if check
            opt.Backend='enumeration';opt.Exhaustive=true;
            [a,EA]=w2b_qp3(un,-A.G,A.c,scale,opt);
            [b,EB]=w2b_qp3(un,-A.G(keep,:),A.c(keep),scale,opt);
            enum_u18(k,:)=a.';enum_u16(k,:)=b.';
            da=norm((a-u)./scale,inf);db=norm((b-v)./scale,inf);
            eg1=kkt(EA);eg2=kkt(EB);
            eraw=A.c(ends)+A.G(ends,:)*b; ert=raw_tolerance(A,ends,un,b,Q);
            eNorm=-eraw./(1+abs(A.c(ends))+abs(A.G(ends,:))*abs(b));
            eChanged=norm((b-a)./scale,inf)>Q.ReferenceMatchTol;
            eCandidate=eChanged && any(eNorm>O.NormalizedPhiTol & -eraw>ert);
            eok=EA.Success && EB.Success && max(eg1,eg2)<=O.KKTGateTol && max(da,db)<=Q.ReferenceMatchTol && ...
                eCandidate==candidate && viol(A,a,1:18)<=O.NormalizedPhiTol && viol(A,b,keep)<=O.NormalizedPhiTol && ...
                EB.Objective<=EA.Objective+Q.ObjectiveTol*(1+EA.Objective);
            M(k,28:32)=[eg1,eg2,da,db,eok];
            if ~eok,status{k}='enumeration_disagreement';continue;end
        end
        if ~base
            status{k}='subset_feasibility_kkt_or_nesting';
        elseif candidate
            M(k,33)=1;status{k}='confirmed';
        elseif changed || any(normalized | outside)
            status{k}='threshold_disagreement';
        elseif any(support(ends))
            status{k}='support_only';
        else
            status{k}='no_detected_change';
        end
    catch ME
        status{k}=['error:' ME.identifier ':' ME.message];
    end
end
t=S.t; path=fullfile(here,'runs',[tag '.mat']);temp=[path '.partial.mat'];
save(temp,'M','fields','status','u0','u18','u16','enum_u18','enum_u16','deletion_dual','keep','t','id','-v7');
movefile(temp,path,'f');
fprintf('A2 %d: %d samples, %d eligible, %d enumerated, %d confirmed\n',id,n,sum(M(:,25)==1),sum(M(:,27)==1),sum(M(:,33)==1));
end

function k=kkt(D)
k=max([D.PrimalNormalized,D.StationarityNormalized,D.DualViolationNormalized,D.ComplementarityNormalized]);
end
function value=viol(A,u,ix)
phi=A.c(ix)+A.G(ix,:)*u;
value=max([0;-phi./(1+abs(A.c(ix))+abs(A.G(ix,:))*abs(u))]);
end
function tol=raw_tolerance(A,ix,un,u,Q)
raw=-A.G(ix,:); scaled=bsxfun(@times,raw,Q.Scale.');
r=hypot(hypot(scaled(:,1),scaled(:,2)),scaled(:,3));nz=r>0;
z=(u-un)./Q.Scale; tol=zeros(numel(ix),1);c=A.c(ix);
d=(c(nz)-raw(nz,:)*un)./r(nz); G=bsxfun(@rdivide,scaled(nz,:),r(nz));
tol(nz)=r(nz).*(Q.FeasTol+Q.RelTol*(1+abs(d)+abs(G)*abs(z)))+64*eps*(1+abs(c(nz))+abs(raw(nz,:))*abs(u));
end
