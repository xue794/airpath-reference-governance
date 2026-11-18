"""Generate named copies; only request, selector and event continuation change."""
from pathlib import Path
import json, hashlib
H=Path(__file__).resolve().parent
S=H.parents[1]/'revision_probe_20261106'/'source_snapshot'

def generate(src, dst, replacements):
    text=src.read_text()
    for old,new in replacements:
        assert text.count(old)==1, (src, old, text.count(old))
        text=text.replace(old,new)
    (H/dst).write_text(text)
    return dict(source=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),output=dst,replacements=replacements)

records=[]
records.append(generate(S/'W2D/w2d_control_snapshot.m','ep6_snapshot.m',[
 ('function S=w2d_control_snapshot(', 'function S=ep6_snapshot('),
 ('K=w2d_case_info(mode,C,O);', "K=struct('addNoise',false,'mode',mode);"),
 ('Gd=w2d_guidance_command(t,r,ze,what,mode,C,O);','Gd=ep6_guidance(t,r,ze,what,C,O);'),
 ('Gov=w2d_governor(xf,u0,C.R,method,O,[]);','Gov=ep6_governor(xf,u0,C.R,method,O);'),
 ('Cmd=viA_feedback_command(y,xf,C.R,C.A);','Cmd=ep6_feedback_command(y,xf,C.R,C.A);')]))
records.append(generate(S/'core/rg_constraint_affine.m','ep6_affine.m',[
 ('function C = rg_constraint_affine(', 'function C = ep6_affine('),
 ('[alpha,ri]=rg_alpha_root(V,X,N,R);',"[alpha,ri]=ep6_alpha_root(V,X,N,R,'reference');")]))
records.append(generate(S/'core/viA_feedback_command.m','ep6_feedback_command.m',[
 ('function C = viA_feedback_command(', 'function C = ep6_feedback_command('),
 ('[alpha,ri]=rg_alpha_root(V,Xc,Nc,R);',"[alpha,ri]=ep6_alpha_root(V,Xc,Nc,R,'feedback');")]))
(H/'generated_source_provenance.json').write_text(json.dumps(records,indent=2))
print('Generated 3 unique functions; original endpoint maps remain unchanged.')
