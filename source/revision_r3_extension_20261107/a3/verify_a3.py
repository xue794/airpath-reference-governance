"""Audit the MATLAB snapshots from raw 36-state output and apply frozen comparisons."""
import csv
import json
from itertools import groupby
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation

from prepare import HERE, ROOT, A1, sha

def mat(path):
    return loadmat(path, simplify_cells=True)

def write_rows(name, rows):
    with (HERE/name).open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

class ForceMap:
    def __init__(self,C):
        self.R=R=C["R"]
        theta=R["P"]["theta"]
        self.cl=theta[[15,16,20,21,22]]
        self.cd=theta[[0,1,4,7,9]]
        self.scale=R["rho"]*R["P"]["S"]/2

    def aero(self,a,V):
        q=self.scale*V**2
        return (q*np.polynomial.polynomial.polyval(a,co) for co in (self.cl,self.cd))

    def M(self,a,V,X):
        lift,drag=self.aero(a,V)
        return lift+(drag+X)*np.tan(a)

    def forces(self,V,gam,nu_v,nu_g,nu_c):
        g=self.R["g"]
        m=self.R["P"]["m"]
        return m*(nu_v+g*np.sin(gam)),m*np.hypot(V*nu_g+g*np.cos(gam),V*np.cos(gam)*nu_c)

    def invert(self,V,X,target):
        lo=np.full_like(V,self.R["alphaLo"])
        hi=np.full_like(V,self.R["alphaHi"])
        assert np.all(self.M(lo,V,X)<=target+1e-10) and np.all(self.M(hi,V,X)>=target-1e-10)
        for _ in range(52):
            mid=(lo+hi)/2
            upper=self.M(mid,V,X)>=target
            lo=np.where(upper,lo,mid)
            hi=np.where(upper,mid,hi)
        return (hi+lo)/2

    def angle_h(self,xf):
        V=xf[:,0];X,N=self.forces(V,xf[:,1],*xf[:,3:6].T)
        a=self.invert(V,X,N)
        lo=self.invert(V,X,self.M(self.R["alphaLo"],V,X)+self.R["DeltaN"])
        hi=self.invert(V,X,self.M(self.R["alphaHi"],V,X)-self.R["DeltaN"])
        return np.column_stack((a-lo,hi-a))

def independent(C,t,x,rows,arm):
    """Use NED wind, quaternion rotation, force polynomials and a separate root inversion."""
    R=C["R"]; F=ForceMap(C); A=C["A"]; W=C["MC"]["wind"]
    wind=W["mean"]+W["amp"]*np.sin(t[:,None]*W["omega"]+W["phase"])
    rot=Rotation.from_quat(x[:,[7,8,9,6]])
    va=rot.apply(x[:,3:6])-wind
    vb=rot.inv().apply(va)
    V=np.linalg.norm(va,axis=1)
    gam=np.arctan2(-va[:,2],np.hypot(va[:,0],va[:,1]))
    chi=np.arctan2(va[:,1],va[:,0]); xf=x[:,19:25]
    vf,gf,cf,nuv,nug,nuc=xf.T
    Xf,Nf=F.forces(vf,gf,nuv,nug,nuc)
    ec=np.arctan2(np.sin(chi-cf),np.cos(chi-cf))
    X,N=F.forces(V,gam,nuv-A["kV"]*(V-vf),nug-A["kGamma"]*(gam-gf),nuc-A["kChi"]*ec)
    aref=F.invert(vf,Xf,Nf); afb=F.invert(V,X,N)
    _,df=F.aero(aref,vf); _,dc=F.aero(afb,V)
    kappa=vb[:,0]/V
    xi=np.column_stack((vf-R["Vb"],gf-R["gammaB"],nuv,nug,nuc))
    rd=np.column_stack((C["VgVarying"]*t,50*np.sin(.01*t),-R["h_ref"]-6*np.sin(.008*t)))
    recomputed=dict(pos_error_m=np.linalg.norm(x[:,:3]-rd,axis=1), kappa=kappa,
        ref_T_N=(df+Xf)/np.cos(aref), T_raw_N=(dc+X)/kappa,T_act_N=x[:,35],
        ref_raw_lower_N=Nf-F.M(R["alphaLo"],vf,Xf),ref_raw_upper_N=F.M(R["alphaHi"],vf,Xf)-Nf,
        fb_raw_lower_N=N-F.M(R["alphaLo"],V,X),fb_raw_upper_N=F.M(R["alphaHi"],V,X)-N,
        ref_alpha_rad=aref,fb_alpha_rad=afb,VB_over_rho=np.einsum("ni,ij,nj->n",xi,R["PB"],xi)/R["rhoB"])
    P=C["P7"]
    recomputed["T_target_N"]=np.clip(recomputed["T_raw_N"],P["Tmin"],P["Tmax"])
    recomputed["T_rate_Nps"]=np.clip((recomputed["T_target_N"]-x[:,35])/P["tauT"],-P["TRateMax"],P["TRateMax"])
    outward=((x[:,35]<=P["Tmin"])&(recomputed["T_rate_Nps"]<0))|((x[:,35]>=P["Tmax"])&(recomputed["T_rate_Nps"]>0))
    recomputed["T_rate_Nps"][outward]=0
    recomputed["T_rate_margin_Nps"]=P["TRateMax"]-abs(recomputed["T_rate_Nps"])
    d=x[:,32:35]
    recomputed["surface_margin_rad"]=np.min(np.minimum(d-P["deltaMin"],P["deltaMax"]-d),axis=1)
    if arm=="Am":
        h=F.angle_h(xf)
        ub=np.column_stack((-R["kVp"]*xi[:,0]-R["kVd"]*nuv,
            -R["kGammaP"]*xi[:,1]-R["kGammaD"]*nug,-R["kChiD"]*nuc))
        direction=np.column_stack((xf[:,3:6],ub))
        eps=1e-5
        # Differentiate the moving gaps directly along the backup vector field.
        sigma=(F.angle_h(xf+eps*direction)-F.angle_h(xf-eps*direction))/(2*eps)+R["kN"]*h
        for i,side in enumerate(["lo","hi"]):
            recomputed["Am_h_"+side+"_rad"]=h[:,i]
            recomputed["Am_backup_"+side+"_radps"]=sigma[:,i]
    errors={col:float(np.max(np.abs(arr-rows[col]))) for col,arr in recomputed.items()}
    assert max(errors.values())<1e-7,errors
    actual_alpha=np.arctan2(vb[:,2],vb[:,0])
    physical=dict(airspeed_margin_mps=np.minimum(V-R["Vmin"],R["Vmax"]-V),
        gamma_margin_rad=R["gammaMax"]-np.abs(gam),
        alpha_margin_rad=np.minimum(actual_alpha-R["alphaPhysicalMin"],R["alphaPhysicalMax"]-actual_alpha),
        T_act_margin_N=np.minimum(x[:,35]-P["Tmin"],P["Tmax"]-x[:,35]),
        surface_margin_rad=rows["surface_margin_rad"],surface_rate_margin_radps=rows["surface_rate_margin_radps"],
        T_rate_margin_Nps=rows["T_rate_margin_Nps"])
    return errors,physical

def clip_spans(t,raw,P):
    categories=np.where(raw<P["Tmin"],-1,np.where(raw>P["Tmax"],1,0))
    spans=[]
    for typ,group in groupby(enumerate(categories),key=lambda z:int(z[1])):
        ix=[i for i,_ in group]
        if typ:
            spans.append(dict(side="upper" if typ>0 else "lower",first_s=float(t[ix[0]]),
                last_s=float(t[ix[-1]]),samples=len(ix)))
    return spans

def read_run(id,arm,level,B):
    base=HERE/"runs"/f"id_{id:03}_{arm}_{level}"
    p=lambda end: base.with_name(base.name+end)
    meta=json.loads(p("_status.json").read_bytes())
    state=mat(p("_states.mat")); diag=mat(p("_samples.mat"))
    t=np.atleast_1d(state["t"]); x=np.atleast_2d(state["x"])
    rows=dict(zip(diag["names"],np.atleast_2d(diag["data"]).T,strict=True))
    assert meta["id"]==id and meta["arm"]==arm and meta["max_step"]=={"fine":.005,"coarse":.01}[level]
    assert state["step"]==meta["max_step"] and meta["status"]==state["status"]
    assert len(t)==meta["n"] and t[-1]==meta["last_saved"]
    assert x.shape==(len(t),36) and np.isfinite(x).all() and np.array_equal(x[0],B["x"][0])
    assert np.array_equal(t,B["t"][:len(t)]) and np.array_equal(t,rows["time_s"])
    assert all(np.isfinite(v).all() for k,v in rows.items() if not k.startswith("Am_") or arm=="Am")
    if meta["status"] in ("completed","archived_completed"): assert t[-1]==120 and len(t)==2401
    if meta["archived_replay"]:
        assert arm=="N10" and level=="coarse" and np.array_equal(B["x"],x)
        assert np.array_equal(B["Y"]["Tcmd"],rows["T_raw_N"])
        assert np.array_equal(B["Y"]["uGov"],np.column_stack([rows[f"u{i}"] for i in [1,2,3]]))
        old=mat(A1/"runs"/f"id_{id:03}_rows.mat")
        for j,side in enumerate(["lower","upper"]):
            for new,name in [("ref","gapRef"),("fb","gapFeedback")]:
                assert np.array_equal(rows[new+"_raw_"+side+"_N"],old[name][:,j])
    errors,physical=independent(B["C"],t,x,rows,arm)
    assert rows["kkt"].max()<=5e-8 and rows["violation"].max()<=1e-8
    return rows,meta,errors,physical

def audit():
    for key,digest in json.loads((HERE/"input_hashes.json").read_bytes()).items():
        assert sha(ROOT/key)==digest,key
    protocol=json.loads((HERE/"protocol.json").read_bytes())
    runs=[]; contrasts=[]; clips=[]; checks=[]; events=[]
    for id in protocol["ids"]:
        B=mat(A1/"runs"/f"id_{id:03}.mat")["S"]; P=B["C"]["P7"]
        cases={(arm,level):read_run(id,arm,level,B) for arm in ["N10","D16","Am"] for level in ["coarse","fine"]}
        last=min(len(r["time_s"]) for r,meta,e,phy in cases.values());start=round(protocol["start_s"][str(id)]/.05)
        minimum={}
        for (arm,level),(s,m,errors,physical) in cases.items():
            t=s["time_s"];sl=slice(0,last)
            primary=s["fb_raw_upper_N"][start:last]
            value=float(np.min(primary)) if primary.size else None
            minimum[arm,level]=value
            k=start+int(np.argmin(primary)) if primary.size else None
            spans=clip_spans(t,s["T_raw_N"],P)
            clips.extend(dict(id=id,arm=arm,step=level,**span) for span in spans)
            row=dict(id=id,arm=arm,step=level,source=m["status"],n=m["n"],original_end_s=t[-1],common_end_s=t[last-1],
                primary_start_s=protocol["start_s"][str(id)],primary_n=len(primary),primary_min_upper_N=value,
                primary_argmin_s=t[k] if k is not None else None,full_common_RMS_m=float(np.sqrt(np.mean(s["pos_error_m"][sl]**2))),
                common_post_RMS_m=float(np.sqrt(np.mean(s["pos_error_m"][start:last]**2))) if primary.size else None,
                clip_samples_common=int(np.count_nonzero(s["T_raw_N"][sl]!=s["T_target_N"][sl])),
                clip_samples_original=sum(span["samples"] for span in spans),last_rhs_s=m["last_rhs"],error_id=m["error_id"])
            for field in ["ref_raw_lower_N","ref_raw_upper_N","fb_raw_lower_N","fb_raw_upper_N","ref_T_N","T_raw_N","T_act_N","T_rate_Nps","VB_over_rho","surface_margin_rad","surface_rate_margin_radps","T_rate_margin_Nps","kkt","violation","backup_violation","Am_backup_lo_radps","Am_backup_hi_radps","Am_h_lo_rad","Am_h_hi_rad"]:
                use=not field.startswith("Am_") or arm=="Am"
                # Ranges describe the complete common interval, including t=0.
                row["min_"+field]=float(min(s[field][sl])) if use else None
                row["max_"+field]=float(max(s[field][sl])) if use else None
            runs.append(row)
            checks.append(dict(id=id,arm=arm,step=level,errors=errors))
            for field,trace in physical.items():
                # Band around numerical zero for counting sampled contacts, in each field's own unit.
                tol=1e-8
                events.append(dict(id=id,arm=arm,step=level,quantity=field,min_margin=min(trace),
                    min_at_s=t[np.argmin(trace)],contact_or_outside_samples=int(np.count_nonzero(trace<=tol)),
                    outside_samples=int(np.count_nonzero(trace < -tol)),contact_tol=tol,
                    status=m["status"],terminal_error=m["error_id"]))
        for other in ["D16","Am"]:
            vals=[minimum[k,v] for k in ["N10",other] for v in ["coarse","fine"]]
            rec=dict(id=id,contrast="N10-"+other,start_s=protocol["start_s"][str(id)],common_end_s=B["t"][last-1])
            if None in vals:
                dc=df=dn=do=gate=None; outcome="window_unavailable"
            else:
                nc,nf,oc,of=vals; dc=nc-oc; df=nf-of; dn=abs(nc-nf); do=abs(oc-of); gate=5*max(1e-7,dn,do)
                passes=dc*df>0 and min(abs(dc),abs(df))>gate
                outcome=("higher" if df>0 else "lower") if passes else "unresolved"
            contrasts.append(dict(**rec,coarse_delta_N=dc,fine_delta_N=df,N10_drift_N=dn,other_drift_N=do,
                five_x_gate_N=gate,classification=outcome))
    write_rows("run_summary.csv",runs);write_rows("contrasts.csv",contrasts);write_rows("model_limits.csv",events)
    with (HERE/"thrust_clipping.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames="id arm step side first_s last_s samples".split());writer.writeheader();writer.writerows(clips)
    summary=dict(cases=len(runs),points=sum(r["n"] for r in runs),completed=sum(bool(r["original_end_s"]==120) for r in runs),
        cases_checked=checks,contrasts=contrasts)
    (HERE/"a3_verified.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
    return runs,contrasts

if __name__=="__main__":
    runs,contrasts=audit()
    print(len(runs),"configurations verified")
    for row in contrasts: print(row)
