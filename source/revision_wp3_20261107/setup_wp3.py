"""Freeze the WP3 postprocessing protocol before reading new numerical results."""
from pathlib import Path
import json, hashlib, shutil, datetime
ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / 'revision_stage1_20261107'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p, v): p.write_text(json.dumps(v, ensure_ascii=False, indent=2), encoding='utf-8')
for n in ['inputs','results','checks','figures','qa','sources']:(ROOT/n).mkdir(exist_ok=True)
catalog=json.loads((BASE/'evidence/trajectory_catalog.json').read_text(encoding='utf-8'))
records=[]
for row in catalog:
    if not row['path']: continue
    assert sha(row['path']) == row['sha256'], row['path']
    records.append(row)
sources=[]
folders=[BASE/'baseline/code/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel',BASE/'baseline/code/P1_Aerospace_W2C_R2018a_v1_0']
for folder in folders:
    for p in sorted(folder.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE/'baseline/code'); dst=ROOT/'sources'/rel
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
        sources.append(dict(source=str(p),snapshot=str(dst),sha256=sha(p)))
extra=BASE/'baseline/transparent_support'
for p in extra.glob('*.m'):
    dst=ROOT/'sources/transparent'/p.name;dst.parent.mkdir(exist_ok=True);shutil.copy2(p,dst)
    sources.append(dict(source=str(p),snapshot=str(dst),sha256=sha(p)))
dump(ROOT/'inputs/trajectories.json',records)
dump(ROOT/'inputs/source_manifest.json',sources)
protocol=dict(created_at=datetime.datetime.now().astimezone().isoformat(),task='WP3',integrations=0,
    cases=['hold','varying','noisy'],records=7,expected_samples=14407,
    missing_records=['Transparent_hold','Transparent_noisy'],missing_policy='WP4 conditional; not filled in WP3',
    model='frozen nominal equations and original gains; R.rhoB=0.10',
    error_box=dict(speed_scale_mps=38.1,coordinates=['eV','38.1*egamma','38.1*wrap(echi)'],norm='infinity',unit='m/s'),
    attitude=dict(primary_deg=0.70,source='original C.tol.tube, conditional parameter, not a proven error bound',
        primary_error='thetaF+thetaR, as in submitted S4',secondary_error='direct SO(3) distance Cc to C; separately reported, justified by triangle inequality'),
    interval=dict(precision_decimal=40,search_cap_mps=1.0,bisections=20,strict_positive=True,
        fixed_constants='original binary64 nominal constants, composite K enclosed from factors',
        geometry='aGamma > aGammaMin and two bank affine margins > 0',
        endpoints='two raw endpoint margins > 0; positive M_alpha over original [0,9deg]',
        joint='geometry + endpoints + kappa>0 + raw pre-actuator thrust strictly in [0,60] N',
        search='three nested sufficient half-widths E, G, J; each final positive value rechecked',
        cap_interpretation='lower sufficient half-width at imposed ceiling; not a maximum'),
    extra_diagnostic='one interval check at each observed error half-width and observed direct attitude distance; posthoc only, no new search',
    tolerance=dict(sample_margin=1e-8,reconstruction_force_N=1e-7,air_state=1e-10,interval_pass=0),
    reports=['all 14407 samples','per-case common sampled-point half-width','positive certificates and inclusion separate','failure attribution','time segments on saved grid only'],
    limits=['No whole-KB uniform transfer certificate','No between-sample guarantee','No error-box invariance or new feedback design','No physical-aircraft safety claim','No mission pass threshold'],
    input_hashes=[dict(path=r['path'],sha256=r['sha256']) for r in records])
target=ROOT/'inputs/protocol.json'
if target.exists():
    old=json.loads(target.read_text(encoding='utf-8'));protocol['created_at']=old['created_at'];assert old==protocol
else:dump(target,protocol)
dump(ROOT/'checks/input_integrity.json',dict(trajectories=7,source_files=len(sources),pass_all=True,protocol_sha256=sha(target)))
print(json.dumps(dict(trajectories=7,sources=len(sources),protocol_sha256=sha(target))))
