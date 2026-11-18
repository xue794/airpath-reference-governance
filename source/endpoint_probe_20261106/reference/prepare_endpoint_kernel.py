from pathlib import Path
import hashlib,json
here=Path(__file__).resolve().parent
source=here.parents[1]/'revision_probe_20261106/source_snapshot/core/rg_constraint_affine.m'
s=source.read_text()
old='''[alpha,ri]=rg_alpha_root(V,X,N,R);
if ~ri.bracket
    error('rg_constraint_affine:root','AoA inversion bracket lost.');
end'''
assert s.count(old)==1
new='''[alpha,ri]=rg_alpha_root(V,X,N,R);
% Numerical extension for locating the original endpoint event only.
% Original endpoint rows above remain unchanged. No post-event trajectory
% is retained or interpreted as valid under the original inversion domain.
if ~ri.bracket
    Re=R;Re.alphaLo=-2*pi/180;Re.alphaHi=12*pi/180;
    [alpha,rie]=rg_alpha_root(V,X,N,Re);
    if ~rie.bracket,error('EndpointProbe:ExtensionRoot','Numerical extension exhausted.');end
end'''
s=s.replace('function C = rg_constraint_affine(x,R)','function C = endpoint_affine_extended(x,R)',1).replace(old,new)
(here/'endpoint_affine_extended.m').write_text(s)
(here/'extended_kernel_provenance.json').write_text(json.dumps({'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'changes':['function rename','root evaluation fallback for numerical event location; original endpoint rows unchanged']},indent=2))
