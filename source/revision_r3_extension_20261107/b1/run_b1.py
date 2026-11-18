"""Run and publish the B1 numerical evidence using the frozen local records."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from interval_check import check
from reserve import uniform,display
from witness import event,old_parameters,legacy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]; A1=HERE.parent/"a1"; A2=HERE.parent/"a2"; A3=HERE.parent/"a3"
WP=ROOT/"revision_wp3_20261107"; CP=ROOT/"contribution_probe_20261106"
CORE=ROOT/"FWUAV/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel"
def digest(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def dump(name,r):
    (HERE/name).write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")
def readj(p):return json.loads(p.read_bytes())
def table(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def write_table(name,rows):
    with (HERE/name).open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def inputs():
    manifests=[A1/"release_manifest.json",A2/"release_manifest.json",A3/"release_manifest.json",
        WP/"artifact_manifest.json",CP/"artifact_manifest.json"]
    known={}
    for m in manifests:
        data=readj(m)
        if isinstance(data,dict):known.update({ROOT/k:v["sha256"] for k,v in data.items()})
        else:known.update({m.parent/row["path"]:row["sha256"] for row in data})
    ids=[int(v["id"]) for v in table(A2/"trial_summary.csv")]
    assert len(set(ids))==len(ids)==17
    sources=["rg_N_map","rg_nominal_aero","rg_params","rg_angle_error","rg_backup_state","viA_feedback_command","rg_constraint_affine"]
    files=[ROOT/"Aerospace_R3端点约束与反馈预留验证方案_v1_0.docx",*[CORE/(s+".m") for s in sources],
        A1/"a1_verified.json",A2/"trial_summary.csv",A2/"event_comparison.csv",A3/"a3_verified.json",
        WP/"inputs/protocol.json",WP/"interval_wp3.py",WP/"results/MatchedQP_hold_export.mat",WP/"results/MatchedQP_hold_transfer.csv",
        CP/"backup/original_certificate/interval_certificate.py",CP/"backup/original_certificate/backup_rho_0p10_certificate.json"]
    for id in ids:files.extend(A1/"runs"/(f"id_{id:03}"+tail) for tail in [".mat","_rows.mat"])
    vendor=sorted(p for p in (WP/"vendor/mpmath").rglob("*.py") if "__pycache__" not in p.parts)
    files+=vendor+[WP/"vendor/mpmath-1.3.0.dist-info/LICENSE"]
    hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}
    for p in files:
        assert p in known and hashes[p.relative_to(ROOT).as_posix()]==known[p],p
    for m in manifests:hashes[m.relative_to(ROOT).as_posix()]=digest(m)
    dump("input_hashes.json",hashes)
    return ids

def kernel(C):
    r=C["R"]; p=r["P"]
    def hx(x):return float(x).hex()
    return dict(Vb=hx(r["Vb"]),gammaB=hx(r["gammaB"]),mass=hx(p["m"]),S=hx(p["S"]),density=hx(r["rho"]),
        g=hx(r["g"]),alphaLo=hx(r["alphaLo"]),alphaHi=hx(r["alphaHi"]),rhoB=hx(r["rhoB"]),DeltaN=hx(r["DeltaN"]),
        kV=hx(C["A"]["kV"]),kGamma=hx(C["A"]["kGamma"]),kChi=hx(C["A"]["kChi"]),
        PB=[[hx(v) for v in row] for row in r["PB"]],
        CL=[hx(p["theta"][i]) for i in [15,16,20,21,22]],CD=[hx(p["theta"][i]) for i in [0,1,4,7,9]])

def replay(C,xf,y,ref,fc):
    R=C["R"];P=R["P"];A=C["A"];m=P["m"];g=R["g"]
    vf,gf,cf,nv,ng,nc=xf.T;v,ga,chi=y.T
    ev=v-vf;eg=ga-gf;ec=np.arctan2(np.sin(chi-cf),np.cos(chi-cf))
    ah=R["alphaHi"];ends=np.array([R["alphaLo"],ah])
    coeff=lambda idx: np.polynomial.polynomial.polyval(ends,P["theta"][idx])
    k=.5*R["rho"]*P["S"];bk=coeff([15,16,20,21,22])+coeff([0,1,4,7,9])*np.tan(ends)
    ref_a=vf*ng+g*np.cos(gf);ref_b=vf*np.cos(gf)*nc;nf=np.hypot(ref_a,ref_b)
    fa=v*(ng-A["kGamma"]*eg)+g*np.cos(ga);fb=v*np.cos(ga)*(nc-A["kChi"]*ec)
    dcos=-2*np.sin(gf+eg/2)*np.sin(eg/2);dsin=2*np.cos(gf+eg/2)*np.sin(eg/2)
    da=ev*ng-v*A["kGamma"]*eg+g*dcos
    db=(ev*np.cos(ga)+vf*dcos)*nc-v*np.cos(ga)*A["kChi"]*ec
    dn=m*(da*(fa+ref_a)+db*(fb+ref_b))/(nf+np.hypot(fa,fb))
    dx=m*(-A["kV"]*ev+g*dsin)
    dm=k*(ev*(2*vf+ev))[:,None]*bk+dx[:,None]*np.tan(ends)
    losses=np.column_stack((dm[:,0]-dn,dn-dm[:,1]))
    def gaps(v,x,n):
        M=k*(v**2)[:,None]*bk+x[:,None]*np.tan(ends)
        return np.column_stack((n-M[:,0],M[:,1]-n))
    rr=gaps(vf,m*(nv+g*np.sin(gf)),m*nf)
    rc=gaps(v,m*(nv-A["kV"]*ev+g*np.sin(ga)),m*np.hypot(fa,fb))
    def err(a):return float(np.max(np.abs(a)))
    errors=[err(rr-ref),err(rc-fc),err(losses-(ref-fc)),err(rr-losses-rc)]
    assert max(errors)<1e-8
    return losses,dict(n=len(vf),ref_replay_N=errors[0],fb_replay_N=errors[1],loss_vs_archive_N=errors[2],
        identity_error_N=errors[3],min_lower_loss_N=float(min(losses[:,0])),max_lower_loss_N=float(max(losses[:,0])),
        min_upper_loss_N=float(min(losses[:,1])),max_upper_loss_N=float(max(losses[:,1])))

def run():
    ids=inputs(); assert readj(A3/"a3_verified.json")["completed"]==18
    mc=[]; nominal=None; witness=None
    first_event=min((row for row in table(A2/"event_comparison.csv") if row["status"]=="confirmed"),
        key=lambda row:(int(row["id"]),float(row["time_s"])))
    for id in ids:
        tag=f"id_{id:03}"; S=loadmat(A1/"runs"/(tag+".mat"),simplify_cells=True)["S"]
        raw=loadmat(A1/"runs"/(tag+"_rows.mat"),simplify_cells=True)
        C=S["C"]; Y=S["Y"]; t=S["t"]
        if nominal is None:nominal=C
        assert kernel(C)==kernel(nominal) and np.array_equal(S["x"][:,19:25],Y["xf"])
        _,tests=replay(C,Y["xf"],np.column_stack([Y[col] for col in ["V","gamma","chi"]]),raw["gapRef"],raw["gapFeedback"])
        mc.append(dict(source="A1",case=tag,**tests))
        if id==int(first_event["id"]):
            idx=int(first_event["index"])-1
            assert id==22 and t[idx]==float(first_event["time_s"])==95.55
            witness=event(S,idx)
    hold=loadmat(WP/"results/MatchedQP_hold_export.mat",simplify_cells=True)
    assert kernel(hold["C"])==kernel(nominal)
    cols={key:k for k,key in enumerate(hold["names"])}
    def take(keys):return hold["D"][:,[cols[k] for k in keys]]
    los,te=replay(hold["C"],take("Vf gammaf chif nuV nuGamma nuChi".split()),take("actual_V actual_gamma actual_chi".split()),
                take(["ref_lower","ref_upper"]),take(["fb_lower","fb_upper"]))
    wp=table(WP/"results/MatchedQP_hold_transfer.csv")
    assert np.max(abs(los-np.array([[float(r["lower_loss_N"]),float(r["upper_loss_N"])] for r in wp])))<1e-8
    mc.append(dict(source="WP3",case="MatchedQP_hold",**te))
    assert sum(v["n"] for v in mc[:-1])==40817 and te["n"]==1201
    cert=display(uniform(nominal)); interval=check(nominal,cert);combined=check(nominal,cert,legacy)
    old=readj(CP/"backup/original_certificate/backup_rho_0p10_certificate.json")
    assert old["certified"] and old["rho"]=="0.1" and witness is not None
    summary=dict(binary_model=kernel(nominal),analytic=cert,independent_binary_model=interval,
        joint_legacy_runtime_intervals=combined,old_parameter_probe=old_parameters(nominal),witness=witness,
        legacy_reserved_ref_lower_N=old["lower_bounds_exact"]["hNminus"],
        legacy_reserved_ref_upper_N=old["lower_bounds_exact"]["hNplus"],MC_saved_points=40817,WP3_saved_points=1201,
        maximum_archived_loss_error_N=max(r["loss_vs_archive_N"] for r in mc),
        maximum_endpoint_replay_error_N=max(max(r["ref_replay_N"],r["fb_replay_N"]) for r in mc),
        b2_status="not run",subdomain_used=False)
    dump("b1_verified.json",summary);write_table("point_checks.csv",mc)
    return summary

if __name__=="__main__":
    outcome=run()
    print("B1 PASS",outcome["analytic"]["radius"],outcome["witness"]["upper_loss_N"])
