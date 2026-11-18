"""Verify the published evidence against the archived samples and source fingerprints."""
import csv,gzip,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from coverage import ROOT,H,A1,A2,A3,WP3,WP4,CP,EP,records,read_data
from boxes import MODEL
from search import writej

def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def source_files():
    b1=H.parent/"b1";r3=ROOT/"revision_r3_value_20261107"
    catalog={}
    manifests=[b1/"release_manifest.json",A3/"release_manifest.json",WP3/"artifact_manifest.json",
               WP4/"artifact_manifest.json",EP/"artifact_manifest.json",CP/"artifact_manifest.json",r3/"release_manifest.json"]
    for m in manifests:
        data=json.loads(m.read_bytes())
        if isinstance(data,list):
            catalog.update({m.parent/r["path"]:r["sha256"] for r in data})
        else:
            catalog.update({ROOT/k:v["sha256"] if isinstance(v,dict) else v for k,v in data.items()})
    # Include B1's verifiable frozen premises, old/new joint interval check and its input manifests.
    files={ROOT/k for k in json.loads((b1/"release_manifest.json").read_bytes())}
    for _,_,paths in records():files.update(paths)
    files.update([r3/"replay/replay_validation.json",r3/"curve/summary.json",
        A3/"run_summary.csv",A3/"contrasts.csv",A3/"model_limits.csv",A3/"thrust_clipping.csv",
        A2/"event_comparison.csv",WP4/"results/summary_table.csv"])
    out={}
    for p in files:
        h=sha(p);assert h==catalog[p],p
        out[p.relative_to(ROOT).as_posix()]=h
    for p in manifests:out[p.relative_to(ROOT).as_posix()]=sha(p)
    out=dict(sorted(out.items()))
    writej(H/"input_hashes.json",out)
    return out

def verify():
    inputs=source_files();summary=json.loads((H/"b2_verified.json").read_bytes());eps=F(summary["radius"])
    frac=lambda v:F.from_float(float(v))
    b=MODEL["binary_model"];pb=[[frac(float.fromhex(x)) for x in row] for row in b["PB"]]
    v0=frac(float.fromhex(b["Vb"]));rho=frac(float.fromhex(b["rhoB"]))
    source_rows=0;near=0;exact_events=0;maxdups=0.
    for group,tag,paths in records():
        C,t,x,y,refA,theta,kappa,raw,gaps,end=read_data(group,tag,paths)
        with gzip.open(H/"coverage"/f"{group}_{tag}.csv.gz","rt") as f:r=list(csv.DictReader(f))
        assert len(r)==len(t);source_rows+=len(r)
        d=np.array([int(s["domain"]) for s in r]);E=np.array([int(s["E"]) for s in r],bool)
        # Direct dense quadratic (not the completed-square interval expression).
        xf=x[:,19:25];xi=xf[:,[0,1,3,4,5]].copy();xi[:,0]-=float(v0)
        ratio=np.einsum("ni,ij,nj->n",xi,C["R"]["PB"],xi)/float(rho)
        nonborder=abs(ratio-1)>1e-11
        assert np.array_equal(d[nonborder],np.where(ratio[nonborder]<1,1,-1))
        select=set(np.flatnonzero(~nonborder).tolist())
        select.update([0,len(t)//2,len(t)-1])
        if group=="A1":
            with (H/"event_coverage.csv").open(encoding="utf-8-sig",newline="") as f:
                idx=[int(s["index"])-1 for s in csv.DictReader(f) if int(s["id"])==int(tag[3:])]
            select.update(idx);exact_events+=len(idx)
        for k in select:
            a,g,_,nv,ng,nc=map(frac,xf[k]);z=[a-v0,g,nv,ng,nc]
            value=sum(pb[i][j]*z[i]*z[j] for i,j in [(0,0),(1,1),(2,2),(3,3),(4,4),(0,2),(2,0),(1,3),(3,1)])
            if d[k]!=0: assert (value<=rho)==(d[k]==1)
        near+=len(select)
        diff=y-xf[:,:3];diff[:,2]=np.arctan2(np.sin(diff[:,2]),np.cos(diff[:,2]))
        widths=np.max(abs(diff)*[1,38.1,38.1],axis=1)
        far=abs(widths-float(eps))>1e-10
        assert np.array_equal(E[far],(d[far]==1)&(widths[far]<float(eps)))
        geom=np.array([bool(int(s["G"])) for s in r]);joint=np.array([bool(int(s["J"])) for s in r])
        assert np.all(~joint|geom) and np.all(~geom|E) and np.all(gaps[E]>0)
        assert np.all(theta[joint]<=math.radians(.7)) and np.all(kappa[joint]>0) and np.all((raw[joint]>0)&(raw[joint]<60))
        assert np.all(np.array([float(s["S4_Tlo_N"]) for s in r])[joint]>0) and np.all(np.array([float(s["S4_Thi_N"]) for s in r])[joint]<60)
        if group=="A3" and tag.endswith("N10_coarse"):
            with gzip.open(H/"coverage"/f"A1_{tag[:6]}.csv.gz","rt") as f:prior=list(csv.DictReader(f))
            for a,p in zip(r,prior,strict=True):
                assert all(a[col]==p[col] for col in ["domain","error_box","E","G","J"])
            maxdups=max(maxdups,max(abs(float(a["theta_sum_deg"])-float(p["theta_sum_deg"])) for a,p in zip(r,prior)))
    assert exact_events==21 and maxdups<1e-5
    writej(H/"sample_review.json",dict(original_input_files=len(inputs),configuration_saved_points=source_rows,
        exact_rational_domain_spot_or_boundary_checks=near,A2_event_indices_checked=exact_events,
        duplicate_A1_A3_attitude_max_difference_deg=maxdups))
    print("Verified sample audit",source_rows,near)
if __name__=="__main__":verify()
