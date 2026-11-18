"""Freeze all third-round scientific artifacts with a verified ZIP manifest."""
from pathlib import Path
import hashlib,json,zipfile

root=Path(__file__).resolve().parent
a=json.loads((root/'comparison_audit.json').read_text())
assert a['complete'] and a['all_available_run_checks_pass']
assert len(a['run_checks'])==6
assert all(x['all_unchanged'] for x in json.loads((root/'prior_round_integrity.json').read_text()))
for name in ['RESULTS_CN.txt','backup/BACKUP_RESULTS_CN.txt','feedback/feedback_transfer_report.txt',
             'aoa/AOA_BASELINE_REVIEW.txt','feedback/posthoc_attitude_summary.json']:
    assert (root/name).is_file(),name
manifest_path=root/'artifact_manifest.json'
files=sorted(p for p in root.rglob('*') if p.is_file() and p!=manifest_path and '__pycache__' not in p.parts)
manifest=[{'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,
           'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
archive=root.with_suffix('.zip')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files+[manifest_path]:z.write(p,arcname=f'{root.name}/{p.relative_to(root).as_posix()}')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for item in manifest:
        assert hashlib.sha256(z.read(f"{root.name}/{item['path']}")).hexdigest()==item['sha256']
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'hashed_artifacts':len(manifest),
                  'integrity':'passed','sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}))
