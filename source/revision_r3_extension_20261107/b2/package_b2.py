"""Seal the portable B2 results together with their frozen numerical inputs."""
import hashlib,json
from pathlib import Path
from zipfile import ZIP_DEFLATED,ZipFile
from verify_saved import source_files,sha

H=Path(__file__).resolve().parent;ROOT=H.parents[1]
def run():
    frozen=source_files();files={ROOT/k for k in frozen}
    for pattern in ["*.py","*.html","*.png","*.svg","*.csv","*.txt","trees/*.csv.gz","coverage/*.csv.gz"]:
        files.update(H.glob(pattern))
    for name in ["input_hashes","search","tree_review","counterexamples","b2_verified","sample_review"]:
        files.add(H/(name+".json"))
    manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)}
    m=H/"release_manifest.json"
    m.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf8")
    target=ROOT/"Aerospace_B2区域预留与轨迹覆盖材料包_v1_0.zip"
    tmp=target.with_suffix(".partial.zip")
    with ZipFile(tmp,"w",ZIP_DEFLATED,compresslevel=2) as z:
        for p in sorted(files|{m}):z.write(p,p.relative_to(ROOT).as_posix())
    with ZipFile(tmp) as z:
        assert set(z.namelist())==set(manifest)|{m.relative_to(ROOT).as_posix()}
        for k,r in manifest.items():
            with z.open(k) as f:assert hashlib.file_digest(f,"sha256").hexdigest()==r["sha256"]
    tmp.replace(target)
    print("Packaged and verified",len(files)+1,"files,",target.stat().st_size,"bytes")
if __name__=="__main__":run()
