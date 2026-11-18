function [R,Bank,Gate] = w2f_preflight(sourceRoot,outdir,O)
%W2F_PREFLIGHT Exact source/data/interface and correctness gate.
root=fileparts(mfilename('fullpath'));
E=load(fullfile(root,'data','w2b_expected_source.mat'));
rows=cell(numel(E.Names),1);
for i=1:numel(E.Names)
    name=E.Names{i};f=fullfile(sourceRoot,name);
    if exist(f,'file')~=2,error('W2F:MissingSource','Missing %s',f);end
    [~,fun]=fileparts(name);actual=which(fun);
    if ~local_same(actual,f)
        error('W2F:PathShadow','%s resolves to %s, expected %s.',fun,actual,f);
    end
    a=fileread(f);a=strrep(a,sprintf('\r\n'),sprintf('\n'));a=strrep(a,sprintf('\r'),sprintf('\n'));
    b=E.Text{i};b=strrep(b,sprintf('\r\n'),sprintf('\n'));b=strrep(b,sprintf('\r'),sprintf('\n'));
    rows{i}=struct('File',name,'ResolvedPath',actual,'SourceTextMatch',strcmp(a,b),'ExpectedRawSHA256',E.SHA256{i});
end
T=table(cellfun(@(r)r.File,rows,'UniformOutput',false),cellfun(@(r)r.ResolvedPath,rows,'UniformOutput',false), ...
    cellfun(@(r)r.SourceTextMatch,rows),cellfun(@(r)r.ExpectedRawSHA256,rows,'UniformOutput',false), ...
    'VariableNames',{'File','ResolvedPath','SourceTextMatch','ExpectedRawSHA256'});
w2f_write_csv(fullfile(outdir,'W2F_source_identity.csv'),T);
assert(all(T.SourceTextMatch),'W2F:SourceChanged','Frozen source identity check failed.');

R=rg_params();assert(abs(R.rhoB-.20)<1e-12,'W2F:Baseline','Expected unchanged source default rhoB=0.20.');
R.rhoB=O.RhoB;
assert(norm(R.nuMax-O.Scale,inf)<1e-14,'W2F:Scale','Frozen QP scale mismatch.');
xt=[R.Vb;R.gammaB;0;0;0;0];Ct=rg_constraint_affine(xt,R);
assert(isequal(Ct.names(:),w2b_row_names()) && isequal(size(Ct.G),[18 3]),'W2F:Rows','Expected locked 18 rows.');

D=load(fullfile(root,'data','w2b_frozen_snapshots.mat'));
if iscell(D.CaseNames),caseNames=D.CaseNames(:);else,caseNames=cellstr(D.CaseNames);end
Bank=struct('X',double(D.X),'CC',double(D.CC),'U0',double(D.U0),'Time',double(D.Time(:)), ...
    'SampleIndex',double(D.SampleIndex(:)),'CaseIndex',double(D.CaseIndex(:)), ...
    'SavedLambda',double(D.SavedLambda(:)),'SavedMinPhi',double(D.SavedMinPhi(:)), ...
    'CaseNames',{caseNames},'RhoB',double(D.RhoB));
assert(abs(Bank.RhoB-O.RhoB)<1e-14,'W2F:BankRho','Snapshot bank rhoB mismatch.');
N=size(Bank.X,1);assert(N==425 && size(Bank.X,2)==6 && size(Bank.U0,2)==3,'W2F:Bank','Unexpected snapshot bank size.');
assert(gcd(O.PermutationStride,N)==1,'W2F:Permutation','PermutationStride must be coprime with bank size.');

Cbank=cell(N,1);uB=zeros(3,N);maxU0Err=0;maxScalarUErr=0;maxScalarLamErr=0;maxQpUErr=0;maxPhiViolation=0;paths=cell(N,1);
for i=1:N
    x=Bank.X(i,:).';cc=Bank.CC(i,:).';u0=Bank.U0(i,:).';
    u0r=rg_nominal_control(x,cc,R);maxU0Err=max(maxU0Err,norm(u0r-u0,inf));
    C=rg_constraint_affine(x,R);Cbank{i}=C;uB(:,i)=rg_backup_control(x,R);
    [us,Ds]=w2f_scalar_selector_core(u0,C,uB(:,i),O.Scale);Gs=rg_scalar_governor(x,u0,R);
    maxScalarUErr=max(maxScalarUErr,norm(us-Gs.u,inf));
    maxScalarLamErr=max(maxScalarLamErr,abs(Ds.lambda-Gs.lambda));
    OQ=O;OQ.FeasibleStart=uB(:,i);[uq,Dq]=w2f_qp_selector_core(u0,C,uB(:,i),OQ);
    Gq=rg_matched_qp(x,u0,R,O,C);
    maxQpUErr=max(maxQpUErr,norm(uq-Gq.u,inf));
    maxPhiViolation=max(maxPhiViolation,max([0;-(C.c+C.G*uq)]));
    paths{i}=Dq.Path;
end
assert(maxU0Err<=O.ReferenceMatchTol,'W2F:U0','Frozen u0 mismatch %.3g.',maxU0Err);
assert(maxScalarUErr<=O.CorrectnessTol && maxScalarLamErr<=O.CorrectnessTol,'W2F:Scalar','Scalar core mismatch.');
assert(maxQpUErr<=O.CorrectnessTol,'W2F:QP','QP core mismatch %.3g.',maxQpUErr);
assert(maxPhiViolation<=1e-8,'W2F:Feasibility','QP feasibility residual too large %.3g.',maxPhiViolation);
Bank.C=Cbank;Bank.uB=uB;Bank.QPPath=paths;

% Case composition and QP path mix.
uc=unique(Bank.CaseIndex(:));caseRows=cell(numel(uc),1);
for k=1:numel(uc)
    j=uc(k);ii=find(Bank.CaseIndex==j);name=Bank.CaseNames{j};
    caseRows{k}=struct('CaseIndex',j,'CaseName',name,'Count',numel(ii),'ScalarActive',sum(Bank.SavedLambda(ii)<1-1e-12));
end
TC=table(cellfun(@(r)r.CaseIndex,caseRows),cellfun(@(r)r.CaseName,caseRows,'UniformOutput',false), ...
    cellfun(@(r)r.Count,caseRows),cellfun(@(r)r.ScalarActive,caseRows), ...
    'VariableNames',{'CaseIndex','CaseName','Count','ScalarActive'});
w2f_write_csv(fullfile(outdir,'W2F_benchmark_bank.csv'),TC);

pnames=unique(paths);pr=cell(numel(pnames),1);
for k=1:numel(pnames),pr{k}=struct('QPPath',pnames{k},'Count',sum(strcmp(paths,pnames{k})));end
TP=table(cellfun(@(r)r.QPPath,pr,'UniformOutput',false),cellfun(@(r)r.Count,pr), ...
    'VariableNames',{'QPPath','Count'});
w2f_write_csv(fullfile(outdir,'W2F_qp_path_mix.csv'),TP);
Gate=struct('Pass',true,'N',N,'MaxU0MatchError',maxU0Err,'MaxScalarUError',maxScalarUErr, ...
    'MaxScalarLambdaError',maxScalarLamErr,'MaxQPUError',maxQpUErr,'MaxQPFeasibilityViolation',maxPhiViolation, ...
    'CaseTable',TC,'PathTable',TP);
end

function same=local_same(a,b)
a=strrep(a,char(92),'/');b=strrep(b,char(92),'/');
if ispc,same=strcmpi(a,b);else,same=strcmp(a,b);end
end
