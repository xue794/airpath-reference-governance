function P=w2d_preflight(sourceRoot,outdir,O)
%W2D_PREFLIGHT Exact frozen source identity, package identity and locked target protocol.
root=fileparts(mfilename('fullpath'));E=load(fullfile(root,'data','w2d_expected_source.mat'));Rows=cell(numel(E.Names),1);
for i=1:numel(E.Names)
 name=E.Names{i};f=fullfile(sourceRoot,name);if exist(f,'file')~=2,error('W2D:MissingSource','Missing %s',f);end
 [~,fun]=fileparts(name);actual=which(fun);if ~local_same(actual,f),error('W2D:PathShadow','%s resolves to %s, expected %s.',fun,actual,f);end
 a=local_norm(fileread(f));b=local_norm(E.Text{i});same=strcmp(a,b);Rows{i}=struct('File',name,'ResolvedPath',actual,'SourceTextMatch',same,'ExpectedRawSHA256',E.SHA256{i});
end
writetable(struct2table([Rows{:}]),fullfile(outdir,'W2D_source_identity.csv'));assert(all(cellfun(@(r)r.SourceTextMatch,Rows)),'W2D:SourceChanged','Frozen source differs from W2-C/W2-D basis.');
for f={'rg_matched_qp','w2b_qp3','w2b_options','w2b_row_names','w2d_control_snapshot','w2d_rhs','w2d_outputs','w2d_reference'}
 actual=which(f{1});expected=fullfile(root,[f{1} '.m']);assert(local_same(actual,expected),'W2D:PackageShadow','%s resolves to %s, expected %s.',f{1},actual,expected);
end
R=rg_params();assert(abs(R.rhoB-.20)<1e-12,'W2D:Baseline','Expected untouched source rhoB=0.20.');R.rhoB=O.RhoB;x=[R.Vb;R.gammaB;0;0;0;0];Ccon=rg_constraint_affine(x,R);assert(isequal(Ccon.names(:),w2b_row_names())&&isequal(size(Ccon.G),[18 3]),'W2D:Rows','Locked 18-row interface changed.');
assert(abs(O.RampStart_s-5)<eps&&abs(O.RampDuration_s-.5)<eps&&abs(O.PulseHold_s-.75)<eps&&abs(O.CaseDuration_s-20)<eps,'W2D:Protocol','Target timing changed.');
[s0,I0]=w2d_smooth_pulse(O.RampStart_s,O.RampStart_s,O.RampDuration_s,O.PulseHold_s);
[s1,I1]=w2d_smooth_pulse(O.RampStart_s+O.RampDuration_s,O.RampStart_s,O.RampDuration_s,O.PulseHold_s);
[s2,I2]=w2d_smooth_pulse(O.RampStart_s+O.RampDuration_s+O.PulseHold_s,O.RampStart_s,O.RampDuration_s,O.PulseHold_s);
[s3,I3]=w2d_smooth_pulse(O.RampStart_s+2*O.RampDuration_s+O.PulseHold_s,O.RampStart_s,O.RampDuration_s,O.PulseHold_s);
assert(abs(s0)<1e-14&&abs(I0)<1e-14&&abs(s1-1)<1e-14&&abs(s2-1)<1e-14&&abs(s3)<1e-14&&abs(I3-(O.RampDuration_s+O.PulseHold_s))<1e-14,'W2D:PulseCheck','C2 pulse endpoint check failed.');
C=cl6_params();C.R.rhoB=O.RhoB;TargetRows=cell(numel(O.Modes),1);for k=1:numel(O.Modes),K=w2d_case_info(O.Modes{k},C,O);TargetRows{k}=struct('Mode',O.Modes{k},'ExpectedRowIndex',K.ExpectedRowIndex,'ExpectedRow',K.ExpectedRow,'TargetV_mps',K.TargetAirPath(1),'TargetGamma_deg',rad2deg(K.TargetAirPath(2)),'TargetChi_deg',rad2deg(K.TargetAirPath(3)),'RampStart_s',K.RampStart_s,'RampDuration_s',K.RampDuration_s,'PulseHold_s',K.PulseHold_s,'Duration_s',K.T);end
writetable(struct2table([TargetRows{:}]),fullfile(outdir,'W2D_predeclared_targets.csv'));
P=struct('Pass',true,'SourceRows',{Rows},'RunRhoB',O.RhoB,'ConstraintRows',18,'SourceModified',false,'TargetRows',{TargetRows},'ComparatorCore','W2-B v1.1 / W2-C v1.0 frozen QP core');
end
function s=local_norm(s),s=strrep(s,sprintf('\r\n'),sprintf('\n'));s=strrep(s,sprintf('\r'),sprintf('\n'));end
function tf=local_same(a,b),a=strrep(a,char(92),'/');b=strrep(b,char(92),'/');if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end,end
