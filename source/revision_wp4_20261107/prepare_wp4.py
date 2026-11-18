"""Inventory existing archives and freeze the two missing WP4 run cells."""
from pathlib import Path
import json,hashlib,shutil,zipfile,re,datetime
H=Path(__file__).resolve().parent;W=H.parent;P=W/'revision_wp3_20261107'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
for n in ['inputs','sources','runs','results','checks','figures','qa','analysis','analysis/inputs','analysis/results','analysis/checks']:(H/n).mkdir(parents=True,exist_ok=True)
archives=[];candidates=[];witnesses=[]
for p in sorted(W.rglob('*.zip')):
    if p.parent!=W and 'FWUAV' not in p.relative_to(W).parts:continue
    row=dict(path=str(p),bytes=p.stat().st_size,matching_entries=[],error=None)
    try:
        with zipfile.ZipFile(p) as z:
            for entry in z.infolist():
                name=entry.filename.lower()
                if name.endswith('.mat') and ('transparent' in name or 'baseline' in name):
                    row['matching_entries'].append(entry.filename)
                    if re.search(r'(transparent|baseline).*(hold|noisy)',name):candidates.append(dict(archive=str(p),entry=entry.filename))
                if 'w1c' in name and name.endswith('.m') and entry.file_size<200000:
                    text=z.read(entry).decode('utf-8',errors='replace')
                    if 'BASELINE_020' in text:
                        lines=text.splitlines();selected=[]
                        for i,line in enumerate(lines):
                            if 'BASELINE_020' in line or 'simulate_complete_6dof_case' in line or 'rhoB' in line:
                                selected.append(dict(line=i+1,text='\n'.join(lines[max(0,i-2):i+3])))
                        witnesses.append(dict(archive=str(p),entry=entry.filename,source_sha256=hashlib.sha256(z.read(entry)).hexdigest(),snippets=selected))
    except Exception as e:row['error']=str(e)
    archives.append(row)
loose=[]
for p in sorted((W/'FWUAV').rglob('*.mat')):
    if re.search(r'transparent|baseline',p.name,re.I):loose.append(dict(path=str(p),bytes=p.stat().st_size))
dump(H/'checks/archive_inventory.json',dict(scope='root-level ZIPs and all FWUAV ZIPs; direct central-directory entries; no recursive nested-ZIP interpretation',
    archives=archives,archive_count=len(archives),loose_candidates=loose,hold_noisy_name_matches=candidates,w1c_witnesses=witnesses))
sources=[]
for p in sorted((P/'sources').rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(P/'sources');d=H/'sources'/rel;d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d)
    sources.append(dict(source=str(p),snapshot=str(d),sha256=sha(p)))
support=W/'revision_probe_20261106/transparent/archive_support'
for filename in ['cl6_control_snapshot_transparent_R2018a.m','cl6_rhs_transparent_R2018a.m','cl6_outputs_transparent_minimal_R2018a.m','simulate_complete_6dof_case_transparent_R2018a.m','run_transparent_governor_baseline_R2018a.m']:
    p=support/filename;d=H/'sources/transparent'/filename
    if d.exists():assert sha(p)==sha(d)
    else:shutil.copy2(p,d);sources.append(dict(source=str(p),snapshot=str(d),sha256=sha(p)))
dump(H/'inputs/source_manifest.json',sources)
copyfiles=['interval_wp3.py','check_interval_wp3.py','assemble_wp3.py']
for filename in copyfiles:shutil.copy2(P/filename,H/'analysis'/filename)
shutil.copytree(P/'vendor',H/'analysis/vendor',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
shutil.copy2(P/'inputs/protocol.json',H/'analysis/inputs/protocol.json')
shutil.copy2(P/'inputs/runtime.json',H/'inputs/prior_analysis_runtime.json')
shutil.copy2(W/'major_revision_plan_20261106/bounded_validation_protocol.txt',H/'inputs/approved_bounded_protocol.txt')
prior=json.loads((P/'inputs/trajectories.json').read_text(encoding='utf-8'))
for r in prior:assert sha(r['path'])==r['sha256']
dump(H/'inputs/reused_trajectories.json',prior)
protocol=dict(created_at=datetime.datetime.now().astimezone().isoformat(),phase='WP4 missing-record branch',
    new_cells=[dict(case='hold',method='Transparent',duration_s=60),dict(case='noisy',method='Transparent',duration_s=120)],
    reused_cells=7,planned_total_cells=9,integrations_this_phase=2,large_MC=0,optional_certified_small=False,optional_V_refinement=False,
    optional_reason='This execution fills only the two explicitly requested missing cells. No radius expansion or new scenario is required. V refinement remains conditional under the approved plan.',
    parameters=dict(rhoB=.10,DeltaN_N=10,DeltaT_N=3,alpha_interval_deg=[0,9],gains='unchanged original cl6_params'),
    numerical=dict(solver='ode45 MATLAB R2018a',RelTol=1e-8,AbsTol=1e-10,MaxStep_s=.01,saved_grid_s=.05,
        progress_callback='records requested saved-grid points, emits progress every 10 simulation seconds; does not change derivatives or solver state'),
    termination=dict(policy='Preserve original errors for raw chart, reference root, feedback root and positive kappa. No model extension or root-domain widening.',
        localization='On exception preserve last successfully returned requested-grid point and attempted RHS time/identifier. No claim of exact event time. No automatic scientific retry.',
        raw_margin='Negative design reserve or backup margin alone is diagnostic and does not terminate; original law must remain defined.'),
    input_identity='Complete C, K, I0, 36-state x0, source paths and source hashes frozen by MATLAB preflight before either integration.',
    noise='Original deterministic time function; seed N/A; no random draws',
    metrics=dict(RMS='sqrt(mean(sum((r-rd)^2,2))) on the original 0.05 s grid',late_windows_s={'hold':30,'noisy':60},
        raw_margin_tolerance=1e-8,phi_normalized_tolerance=1e-8,backup_tolerance=1e-8,hold_state_match_tolerance=1e-7,
        transfer='Unchanged WP3 E/G/J interval evaluator; speed scale 38.1 m/s; 0.70 degree S4 attitude sum; cap 1 m/s; 20 bisections'),
    negative_findings='Retain all failures, all partial horizons and all constraint/tracking tradeoffs. No tuning for a positive conclusion.',
    controller_change='Only Transparent reference selector u=u0; original downstream projection, actuator dynamics and inputs retained.')
p=H/'inputs/protocol_manifest.json'
if p.exists():
    old=json.loads(p.read_text(encoding='utf-8'));protocol['created_at']=old['created_at'];assert protocol==old
else:dump(p,protocol)
dump(H/'checks/preparation.json',dict(source_files=len(sources),archives=len(archives),archive_errors=sum(x['error'] is not None for x in archives),
    baseline_name_matches=len(candidates),protocol_sha256=sha(p),source_manifest_sha256=sha(H/'inputs/source_manifest.json')))
print(json.dumps(dict(archives=len(archives),errors=sum(x['error'] is not None for x in archives),candidate_names=candidates,source_files=len(sources)),ensure_ascii=False))
