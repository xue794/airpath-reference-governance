function P=w2e_preflight(sourceRoot,outdir,O)
%W2E_PREFLIGHT Exact source identity, exact saved-sample identity, and protocol lock.
root=fileparts(mfilename('fullpath'));E=load(fullfile(root,'data','w2e_expected_source.mat'));Rows=cell(numel(E.Names),1);
for i=1:numel(E.Names)
 name=E.Names{i};f=fullfile(sourceRoot,name);if exist(f,'file')~=2,error('W2E:MissingSource','Missing frozen source file: %s',f);end;[~,fun]=fileparts(name);actual=which(fun);if ~local_same(actual,f),error('W2E:PathShadow','%s resolves to %s, expected %s.',fun,actual,f);end;a=local_norm(fileread(f));b=local_norm(E.Text{i});same=strcmp(a,b);Rows{i}=struct('File',name,'ResolvedPath',actual,'SourceTextMatch',same,'ExpectedRawSHA256',E.SHA256{i});
end
writetable(struct2table([Rows{:}]),fullfile(outdir,'W2E_source_identity.csv'));assert(all(cellfun(@(r)r.SourceTextMatch,Rows)),'W2E:SourceChanged','Frozen source differs from the W2-E basis.');
for f={'rg_matched_qp','w2b_qp3','w2b_options','w2b_row_names','w2e_governor','w2e_control_snapshot','w2e_rhs','w2e_outputs'}
 actual=which(f{1});expected=fullfile(root,[f{1} '.m']);assert(local_same(actual,expected),'W2E:PackageShadow','%s resolves to %s, expected %s.',f{1},actual,expected);
end
Rmc=rev41_mc_params();assert(Rmc.N==O.ExpectedOldN&&Rmc.seed==O.ExpectedOldSeed&&strcmp(Rmc.mode,O.Mode)&&strcmp(Rmc.initStrategy,'backup_feasible'),'W2E:MCProtocol','Frozen Rev.4.1.1 MC definition changed.');
Cbase=cl6_params();assert(abs(Cbase.R.rhoB-.20)<1e-12,'W2E:Baseline','Untouched source must retain historical rhoB=0.20.');Cbase.R.rhoB=O.RhoB;
x=[Cbase.R.Vb;Cbase.R.gammaB;0;0;0;0];Ccon=rg_constraint_affine(x,Cbase.R);assert(isequal(Ccon.names(:),w2b_row_names())&&isequal(size(Ccon.G),[18 3]),'W2E:Rows','Locked 18-row interface changed.');
% Load the complete saved uncertainty records; do not rely on seed alone.
L=load(O.SampleSource,'Results');assert(isfield(L,'Results')&&isfield(L.Results,'Samples'),'W2E:SampleFile','Frozen sample source does not contain Results.Samples.');Samples=L.Results.Samples;assert(numel(Samples)==O.N,'W2E:SampleCount','Frozen sample count is not 500.');idx=[Samples.index];assert(isequal(idx(:),(1:O.N).'),'W2E:SampleIndex','Frozen sample indices are not exactly 1:500.');
% Independent deterministic regeneration is an identity audit only. Save/restore caller RNG state.
rs=rng;cr=onCleanup(@()rng(rs));Generated=rev41_generate_samples(Rmc);SampleCompare=w2e_compare_samples(Samples,Generated);clear cr;rng(rs);assert(SampleCompare.ExactMatch,'W2E:SampleMismatch','Saved full Samples do not exactly match the frozen generator; do not continue.');
sampleHash=w2e_sha256(O.SampleSource);SampleRow=struct('SampleFile',O.SampleSource,'SHA256',sampleHash,'N',numel(Samples),'Seed',Rmc.seed,'ExactRegenerationMatch',SampleCompare.ExactMatch,'MismatchCount',SampleCompare.MismatchCount);writetable(struct2table(SampleRow),fullfile(outdir,'W2E_sample_identity.csv'));
P=struct('Pass',true,'RunRhoB',O.RhoB,'ConstraintRows',18,'SourceRows',{Rows},'Rmc',Rmc,'Cbase',Cbase,'Samples',Samples,'SampleHash',sampleHash,'SampleCompare',SampleCompare,'SourceModified',false);
end
function s=local_norm(s),s=strrep(s,sprintf('\r\n'),sprintf('\n'));s=strrep(s,sprintf('\r'),sprintf('\n'));end
function tf=local_same(a,b),a=strrep(a,char(92),'/');b=strrep(b,char(92),'/');if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end,end
