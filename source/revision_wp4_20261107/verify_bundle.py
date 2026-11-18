"""Verify packaged relative-path hashes without requiring the original workspace."""
from pathlib import Path
import hashlib,json
H=Path(__file__).resolve().parent
items=json.loads((H/'artifact_manifest.json').read_text(encoding='utf-8'))
for r in items:
    p=H/r['path'];assert p.is_file(),r['path']
    assert p.stat().st_size==r['bytes'],r['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']
print(json.dumps(dict(bundle_files_verified=len(items),all_passed=True)))
