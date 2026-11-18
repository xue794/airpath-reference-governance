"""Create isolated AoA comparison harness from the audited second-round runner."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
old = root.parent / 'endpoint_probe_20261106' / 'full6dof'
dest = root / 'full6dof'
dest.mkdir(exist_ok=True)
manifest = []
for name in ['ep6_affine', 'ep6_alpha_root', 'ep6_feedback_command', 'ep6_guidance',
             'ep6_raw_margins', 'ep6_reference', 'ep6_snapshot_warm', 'ep6_governor_warm',
             'run_ep6_warm']:
    source = old / (name + '.m')
    body = source.read_text(encoding='utf-8').replace('ep6_', 'cp6_').replace('EP6', 'CP6')
    if name == 'run_ep6_warm':
        body = body.replace("methods={'T','N0','N10'}", "methods={'A0','Am'}")
        body = body.replace("addpath(here);", "addpath(here);addpath(fullfile(fileparts(here),'aoa'));")
        body = body.replace("fullfile(outdir,'all_summary.json')", "fullfile(outdir,[strjoin(methods,'_') '_batch_summary.json'])")
    if name == 'ep6_governor_warm':
        body = body.replace("if strcmp(method,'N0'),R.DeltaN=0;end", "if strcmp(method,'N0')||strcmp(method,'A0'),R.DeltaN=0;end")
        before = "if strcmp(method,'T'),ix=[1:13 16:17];else,ix=1:17;end\nA=-C.G(ix,:);b=C.c(ix);Q=O.QP;Q.FeasibleStart=rg_backup_control(x,R);"
        after = """if strcmp(method,'A0')||strcmp(method,'Am')
    ix=1:17;Ac=cp_aoa_rows(x,R,method,C.alpha);
    cc=C.c(ix);GG=C.G(ix,:);cc(14:15)=Ac.c;GG(14:15,:)=Ac.G;
elseif strcmp(method,'T')
    ix=[1:13 16:17];cc=C.c(ix);GG=C.G(ix,:);
else
    ix=1:17;cc=C.c(ix);GG=C.G(ix,:);
end
A=-GG;b=cc;Q=O.QP;Q.FeasibleStart=rg_backup_control(x,R);"""
        assert before in body
        body = body.replace(before, after)
        body = body.replace("'phi',C.c(ix)+C.G(ix,:)*u,'rows',ix", "'phi',cc+GG*u,'rows',ix")
    target = dest / (name.replace('ep6_', 'cp6_') + '.m')
    target.write_text(body, encoding='utf-8')
    manifest.append({'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                     'generated': target.name, 'generated_sha256': hashlib.sha256(target.read_bytes()).hexdigest()})
(dest/'generated_provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps({'generated': len(manifest), 'folder': str(dest)}))
