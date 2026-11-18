"""Package exactly the B1 analysis and its verified dependencies."""
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED,ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def sha(stream):
    return hashlib.file_digest(stream,"sha256").hexdigest()

def package():
    locked=json.loads((HERE/"input_hashes.json").read_bytes())
    files={ROOT/key for key in locked}
    for pattern in ["*.py","*.csv","*.html","*.txt"]:
        files.update(HERE.glob(pattern))
    files.update([HERE/"b1_verified.json", HERE/"input_hashes.json"])
    manifest={}
    for p in sorted(files):
        key=p.relative_to(ROOT).as_posix()
        with p.open("rb") as stream: fingerprint=sha(stream)
        if key in locked: assert locked[key]==fingerprint,key
        manifest[key]={"bytes":p.stat().st_size,"sha256":fingerprint}
    m=HERE/"release_manifest.json"
    m.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    target=ROOT/"Aerospace_B1端点预留与误差半宽材料包_v1_0.zip"
    temp=target.with_suffix(".partial.zip")
    with ZipFile(temp,"w",ZIP_DEFLATED,compresslevel=3) as z:
        for p in sorted(files|{m}):z.write(p,p.relative_to(ROOT).as_posix())
    with ZipFile(temp) as z:
        assert set(z.namelist())==set(manifest)|{m.relative_to(ROOT).as_posix()}
        for key,row in manifest.items():
            with z.open(key) as f:assert sha(f)==row["sha256"],key
    temp.replace(target)
    print("B1 release verified:",len(files)+1,"files,",target.stat().st_size,"bytes")

if __name__=="__main__":
    package()
