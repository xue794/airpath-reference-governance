"""Seal A3 inputs, verify both earlier release manifests and generate the two adapters."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PKG = ROOT / "FWUAV/P1_Aerospace_W2E_R2018a_v1_0"
CORE = ROOT / "FWUAV/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel"
A1, A2 = (HERE.parent / sub for sub in ["a1", "a2"])

def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def prepare():
    for sub in ["runs","identity"]:
        (HERE/sub).mkdir(exist_ok=True)
    for src,dst,subs in [
        ("w2e_rhs","a3_rhs",[("xdot=w2e_rhs(","xdot=a3_rhs("),("S=w2e_control_snapshot(","S=a3_snapshot(")]),
        ("w2e_control_snapshot","a3_snapshot",[("S=w2e_control_snapshot(","S=a3_snapshot("),
          ("Gov=w2e_governor(xf,u0,C.R,method,O,[]);","Gov=a3_governor(xf,u0,C.R,method,O);")])]:
        body=(PKG/(src+".m")).read_text()
        for before,after in subs:
            assert body.count(before)==1
            body=body.replace(before,after)
        (HERE/(dst+".m")).write_text(body,encoding="utf-8",newline="\n")
    a=json.loads((A2/"a2_verified.json").read_bytes())
    with (A2/"trial_summary.csv").open(encoding="utf-8-sig") as f:
        cases=[(int(row["id"]),float(row["first_confirmed_s"])) for row in csv.DictReader(f) if int(row["confirmed"])>0]
    ids=sorted(i for i,t in cases)[:3]
    assert ids==a["a3_ids"]==[22,32,33] and [a["eligible"],a["confirmed"],a["exceptions"]]==[40817,21,0]
    p=dict(ids=ids,start_s={i:dict(cases)[i] for i in ids},grid_s=.05,end_s=120,steps_s=[.01,.005],
        metric="minimum raw upper feedback endpoint gap in N over saved grid from first A2 confirmation inclusive to common end across all six runs per ID",
        contrast_order=["N10-D16","N10-Am"],resolvable="same nonzero sign and min absolute coarse/fine contrast strictly greater than 5*max(abs(coarse-fine) for each of the two arms,1e-7 N)",
        reports="Each original horizon, common interval, lower/upper reference/feedback gaps, raw and actual thrust, actuator limits, backup and RMS. Am witness only on its saved states.")
    (HERE/"protocol.json").write_text(json.dumps(p,indent=2))
    inputs=[*CORE.glob("*.m"),*PKG.glob("*.m"),PKG/"data/w2e_expected_source.mat",
        ROOT/"Aerospace_R3端点约束与反馈预留验证方案_v1_0.docx",
        A1/"a1_verified.json",A2/"a2_verified.json",A2/"trial_summary.csv",A2/"event_comparison.csv"]
    for i in ids:
        tag=f"id_{i:03}"
        inputs.extend([A1/"runs"/(tag+".mat"),A1/"runs"/(tag+"_rows.mat"),A2/"runs"/(tag+".mat")])
    release={}
    for dr in (A1,A2):
        release.update(json.loads((dr/"release_manifest.json").read_bytes()))
    for f in inputs:
        assert sha(f)==release[f.relative_to(ROOT).as_posix()]["sha256"],f
    inputs.extend([*HERE.glob("*.m"),HERE/"protocol.json",HERE/"prepare.py",ROOT/"contribution_probe_20261106/aoa/cp_aoa_rows.m"])
    (HERE/"input_hashes.json").write_text(json.dumps({f.relative_to(ROOT).as_posix():sha(f) for f in sorted(inputs)},indent=2))
    print("Verified A2 and sealed",len(inputs),"inputs.")

if __name__=="__main__":
    prepare()
