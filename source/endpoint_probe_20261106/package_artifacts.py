"""Archive the completed exploratory validation with hashes of every artifact."""
from pathlib import Path
import hashlib
import json
import zipfile

root = Path(__file__).resolve().parent
audit = json.loads((root / 'review/final_validation.json').read_text())
assert audit['all_expected_results_present'] and audit['all_available_run_checks_pass']
integrity = json.loads((root / 'first_round_integrity.json').read_text())
assert integrity['all_unchanged']
manifest_path = root / 'artifact_manifest.json'
files = sorted(p for p in root.rglob('*') if p.is_file()
               and p != manifest_path and '__pycache__' not in p.parts)
manifest = [{'path': p.relative_to(root).as_posix(), 'bytes': p.stat().st_size,
             'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
archive = root.with_suffix('.zip')
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in files + [manifest_path]:
        z.write(p, arcname=f'{root.name}/{p.relative_to(root).as_posix()}')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for item in manifest:
        data = z.read(f"{root.name}/{item['path']}")
        assert hashlib.sha256(data).hexdigest() == item['sha256']
print(json.dumps({'archive': str(archive), 'bytes': archive.stat().st_size,
                  'hashed_artifacts': len(manifest), 'zip_integrity': 'passed',
                  'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()}))
