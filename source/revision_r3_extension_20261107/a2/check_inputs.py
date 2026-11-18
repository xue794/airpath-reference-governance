"""Verify the archived inputs before running the A2 MATLAB replay."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
A1 = HERE.parent / "a1"
CORE = ROOT / "FWUAV/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel"
PKG = ROOT / "FWUAV/P1_Aerospace_W2E_R2018a_v1_0"
ADMISSION = ROOT / "FWUAV/P1_Aerospace_W2B_R2018a_v1_1/w2b_admission.m"

def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()

def check_inputs():
    original = json.loads((A1/"release_manifest.json").read_bytes())
    summary = json.loads((A1/"a1_verified.json").read_bytes())
    assert (summary["accepted"],summary["grid_points"],summary["upper_support_samples"])==(17,40817,21)
    inputs = sorted([*CORE.glob("*.m"), *PKG.glob("*.m"), *list((A1/"runs").glob("*.mat")),
                     PKG/"data/w2e_expected_source.mat", PKG/"data/Rev41Results_FROZEN_SAMPLE_SOURCE.mat",
                     A1/"a1_verified.json", A1/"event_checks.csv", A1/"identity/frozen.mat",
                     ROOT/"Aerospace_R3端点约束与反馈预留验证方案_v1_0.docx"])
    verified={}
    for p in inputs:
        name=p.relative_to(ROOT).as_posix()
        verified[name]=digest(p)
        assert original[name]["sha256"]==verified[name],name
    assert len(list((A1/"runs").glob("*.mat")))==34
    assert digest(ADMISSION)==digest(ROOT/"revision_stage1_20261107/baseline/code/P1_Aerospace_W2B_R2018a_v1_1/w2b_admission.m")
    verified[ADMISSION.relative_to(ROOT).as_posix()]=digest(ADMISSION)
    return verified

if __name__=="__main__":
    items=check_inputs()
    (HERE/"input_hashes.json").write_text(json.dumps(items,ensure_ascii=False,indent=2),encoding="utf-8")
    print("Verified",len(items),"A2 inputs.")
