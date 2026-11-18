"""Refresh release metadata after an intentional edit; protect frozen inputs."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PARTS = {'.git', 'generated', '__pycache__', '.pytest_cache'}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    provenance = json.loads((ROOT/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    for item in provenance['copied_files']:
        path = ROOT/item['file']
        if not path.is_file() or sha256(path) != item['sha256']:
            raise SystemExit('Frozen provenance mismatch; refusing to refresh: ' + item['file'])
    records = []
    for path in sorted(ROOT.rglob('*')):
        relative = path.relative_to(ROOT)
        if (not path.is_file() or path.name.endswith(('.pyc', '.pyo'))
                or any(part in EXCLUDED_PARTS for part in relative.parts)
                or relative.as_posix() == 'RELEASE_MANIFEST.json'):
            continue
        records.append({'file': relative.as_posix(), 'bytes': path.stat().st_size,
                        'sha256': sha256(path)})
    target = ROOT/'RELEASE_MANIFEST.json'
    target.write_text(json.dumps({'version':'v1.0.0', 'files':records}, indent=2) + '\n',
                      encoding='utf-8')
    print('Manifest refreshed:', len(records), 'files; frozen provenance unchanged.')
    print('Run python scripts/verify_release.py before publishing.')


if __name__ == '__main__':
    main()
