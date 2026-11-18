from pathlib import Path
import json,hashlib,shutil
H=Path(__file__).resolve().parent;P=H.parent/'revision_wp3_20261107'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
new=[]
for mode,T in [('hold',60),('noisy',120)]:
    status=json.loads((H/'runs'/f'Transparent_{mode}_status.json').read_text(encoding='utf-8'))
    assert status['Status'] in ['completed','raw_command_chart_exit','reference_inversion_domain_exit','feedback_inversion_domain_exit','thrust_projection_domain_exit','integrator_stopped_early','unclassified_exception']
    path=H/'runs'/f'Transparent_{mode}.mat'
    new.append(dict(case=mode,method='Transparent',duration_s=T,path=str(path),sha256=sha(path),status=status['Status'],new_integration=True))
old=json.loads((H/'inputs/reused_trajectories.json').read_text(encoding='utf-8'))
for r in old:assert sha(r['path'])==r['sha256']
dump(H/'inputs/new_trajectories.json',new);dump(H/'inputs/all_trajectories.json',old+new)
dump(H/'analysis/inputs/trajectories.json',new)
shutil.copy2(H/'inputs/source_manifest.json',H/'analysis/inputs/source_manifest.json')
definition=json.loads((P/'inputs/protocol.json').read_text(encoding='utf-8'))
definition['task']='WP4 application of unchanged WP3 numerical definitions';definition['records']=2;definition['cases']=['hold','noisy'];definition['expected_samples']=3602
definition['reports']=['new two cells only; original seven results separately reused'];definition['input_hashes']=[dict(path=x['path'],sha256=x['sha256']) for x in new]
definition['inherited_definition_sha256']=sha(P/'inputs/protocol.json');definition['missing_records']=[]
dump(H/'analysis/inputs/protocol.json',definition)
# Only wrapper/root and partial-horizon assertion change; numerical exporter is unchanged.
src=(P/'export_wp3.m').read_text(encoding='utf-8');dst=src.replace('function export_wp3()','function export_wp4()').replace("'trajectories.json'","'new_trajectories.json'")
dst=dst.replace("assert(size(x,2)==36 && n==size(x,1) && abs(t(end)-r.duration_s)<1e-10);","assert(size(x,2)==36 && n==size(x,1) && t(end)<=r.duration_s+1e-10);")
dst=dst.replace('WP3_EXPORT_COMPLETE','WP4_EXPORT_COMPLETE');(H/'export_wp4.m').write_text(dst,encoding='utf-8')
check=(P/'check_interval_wp3.py').read_text(encoding='utf-8').replace('input_trajectories_unchanged=7','input_trajectories_unchanged=len(cat)')
(H/'analysis/check_interval_wp3.py').write_text(check,encoding='utf-8')
prior=H/'analysis/prior_results';prior.mkdir(exist_ok=True)
refs=[]
for p in [P/'results/transfer_summary.json',P/'results/joint_failure_attribution.json',P/'checks/independent_verification.json']+list((P/'results').glob('*_transfer.csv')):
    d=prior/p.name;shutil.copy2(p,d);refs.append(dict(source=str(p),copy=str(d),sha256=sha(p)))
dump(H/'inputs/prior_analysis_manifest.json',refs)
derivation=dict(original_exporter_sha256=sha(P/'export_wp3.m'),adapted_exporter_sha256=sha(H/'export_wp4.m'),
    unchanged_interval_evaluator=sha(P/'interval_wp3.py')==sha(H/'analysis/interval_wp3.py'),
    changes=['export function renamed and catalog points to two new files','partial end-time accepted without extending data','independent checker input-count label made dynamic'],
    numerical_thresholds_changed=False)
assert derivation['unchanged_interval_evaluator'];dump(H/'checks/analysis_derivation.json',derivation)
print(json.dumps(dict(new_cells=len(new),reused=len(old),total=len(new)+len(old),all_new_completed=all(x['status']=='completed' for x in new))))
