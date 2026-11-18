function P=w2c_preflight(sourceRoot,outdir,O)
%W2C_PREFLIGHT Exact uploaded-source identity, path resolution, and frozen protocol checks.
root=fileparts(mfilename('fullpath'));E=load(fullfile(root,'data','w2c_expected_source.mat'));Rows=cell(numel(E.Names),1);
for i=1:numel(E.Names)
 name=E.Names{i};f=fullfile(sourceRoot,name);if exist(f,'file')~=2,error('W2C:MissingSource','Missing %s',f);end
 [~,fun]=fileparts(name);actual=which(fun);if ~local_same(actual,f),error('W2C:PathShadow','%s resolves to %s, expected %s.',fun,actual,f);end
 a=fileread(f);a=strrep(a,sprintf('\r\n'),sprintf('\n'));a=strrep(a,sprintf('\r'),sprintf('\n'));
 b=E.Text{i};b=strrep(b,sprintf('\r\n'),sprintf('\n'));b=strrep(b,sprintf('\r'),sprintf('\n'));same=strcmp(a,b);
 Rows{i}=struct('File',name,'ResolvedPath',actual,'SourceTextMatch',same,'ExpectedRawSHA256',E.SHA256{i});
end
local_write_csv(fullfile(outdir,'W2C_source_identity.csv'),Rows);
assert(all(cellfun(@(r)r.SourceTextMatch,Rows)),'W2C:SourceChanged','Frozen source text differs from uploaded W2-C basis.');
for f={'rg_matched_qp','w2b_qp3','w2b_options','w2b_row_names','w2c_control_snapshot','w2c_rhs','w2c_outputs'}
 actual=which(f{1});expected=fullfile(root,[f{1} '.m']);assert(local_same(actual,expected),'W2C:PackageShadow','%s resolves to %s, expected %s.',f{1},actual,expected);
end
R=rg_params();assert(abs(R.rhoB-.20)<1e-12,'W2C:Baseline','Expected untouched source rhoB=0.20.');R.rhoB=O.RhoB;
x=[R.Vb;R.gammaB;0;0;0;0];C=rg_constraint_affine(x,R);assert(isequal(C.names(:),w2b_row_names()) && isequal(size(C.G),[18 3]),'W2C:Rows','Locked 18-row interface changed.');
assert(norm(O.QP.Scale-[2;8*pi/180;12*pi/180],inf)<1e-14,'W2C:Scale','QP metric changed.');
P=struct('Pass',true,'SourceRows',{Rows},'RunRhoB',O.RhoB,'ConstraintRows',18,'SourceModified',false,'ComparatorCore','W2B v1.1 byte-for-byte copy');
end
function tf=local_same(a,b)
a=strrep(a,char(92),'/');b=strrep(b,char(92),'/');if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end
end
function local_write_csv(file,C)
if isempty(C),return;end;T=struct2table([C{:}]);writetable(T,file);
end
