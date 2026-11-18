from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'RELEASE_MANIFEST.json').read_text(encoding='utf-8'))
for record in manifest['files']:
    p=ROOT/record['file']
    assert p.is_file(),str(p)
    assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'],str(p)
print('PASS:',len(manifest['files']),'release file hashes')
