"""Read-only audit of the generated full-6DOF comparison harness."""
from pathlib import Path
import hashlib,json,difflib
root=Path(__file__).resolve().parent.parent
new=root/'full6dof'
manifest=json.loads((new/'generated_provenance.json').read_text())
items=[];diffs=[]
for row in manifest:
    src=Path(row['source']);dst=new/row['generated']
    expected=src.read_text(encoding='utf-8').replace('ep6_','cp6_').replace('EP6','CP6')
    actual=dst.read_text(encoding='utf-8')
    intended_changes=src.stem in {'ep6_governor_warm','run_ep6_warm'}
    unchanged_except_names=expected==actual
    items.append({'file':dst.name,
      'source_hash_matches_manifest':hashlib.sha256(src.read_bytes()).hexdigest()==row['source_sha256'],
      'generated_hash_matches_manifest':hashlib.sha256(dst.read_bytes()).hexdigest()==row['generated_sha256'],
      'only_identifier_names_changed':unchanged_except_names,
      'additional_changes_expected':intended_changes})
    assert unchanged_except_names or intended_changes, f'Unexpected changes in {dst}'
    if not unchanged_except_names:
        diffs.extend(difflib.unified_diff(expected.splitlines(),actual.splitlines(),
                      fromfile=src.name+' (names normalized)',tofile=dst.name,lineterm=''))
assert all(x['source_hash_matches_manifest'] and x['generated_hash_matches_manifest'] for x in items)
out={'all_manifest_hashes_match':True,'files':items,
     'aoa_kernel_sha256':hashlib.sha256((root/'aoa/cp_aoa_rows.m').read_bytes()).hexdigest()}
(root/'aoa/harness_source_audit.json').write_text(json.dumps(out,indent=2))
(root/'aoa/harness_intended_diff.txt').write_text('\n'.join(diffs),encoding='utf-8')
print(json.dumps(out,indent=2))
