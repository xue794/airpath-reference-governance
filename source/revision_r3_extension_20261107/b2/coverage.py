"""Classify frozen sampled states, retaining membership ambiguity and S4 denominators."""
from concurrent.futures import ProcessPoolExecutor
import csv,gzip,io,json,math,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from boxes import H,MODEL,Kernel,I,pt,S
from search import writej
ROOT=H.parents[1];WP3=ROOT/"revision_wp3_20261107"; WP4=ROOT/"revision_wp4_20261107"
CP=ROOT/"contribution_probe_20261106"; EP=ROOT/"endpoint_probe_20261106"; A1=H.parent/"a1";A2=H.parent/"a2";A3=H.parent/"a3"
sys.path.insert(0,str(H.parent/"b1"))
from run_b1 import kernel
sys.path.insert(0,str(WP3))
import interval_wp3 as wp

def records():
    jobs=[]
    for p in sorted([*(WP3/"results").glob("*_export.mat"),*(WP4/"results").glob("*_export.mat")]):
        jobs.append(("WP4",p.stem.replace("_export",""),[p]))
    for p in sorted((A1/"runs").glob("id_???.mat")):
        jobs.append(("A1",p.stem,[p,p.with_name(p.stem+"_rows.mat")]))
    for folder in [EP/"full6dof",CP/"full6dof"]:
        for speed in [31,32]:
            name=f"V{speed}_C60_"+("warm_fine" if folder.parent==EP else "fine")
            for p in sorted((folder/name).glob("*.mat")):
                jobs.append(("stress",f"V{speed}_{p.stem}",[p]))
    for p in sorted((A3/"runs").glob("*_states.mat")):
        tag=p.name.replace("_states.mat","")
        jobs.append(("A3",tag,[p,p.with_name(tag+"_samples.mat"),A1/"runs"/(tag[:6]+".mat")]))
    return jobs

def load(p):return loadmat(p,simplify_cells=True)
def ascols(d,cols):return np.column_stack([d[k] for k in cols])
def named(j,k,fields):return dict(zip(j[fields],j[k].T,strict=True))
def air(x,wind):
    body=Rotation.from_quat(x[:,[7,8,9,6]])
    va=body.apply(x[:,3:6])-wind
    return np.column_stack((np.linalg.norm(va,axis=1),np.arctan2(-va[:,2],np.hypot(va[:,0],va[:,1])),
                            np.arctan2(va[:,1],va[:,0]))),body
def attitude(C,x,y,body,ac):
    R,A=C["R"],C["A"];xf=x[:,19:25];v,g,chi=y.T
    err=chi-xf[:,2];ec=np.arctan2(np.sin(err),np.cos(err))
    ag=v*(xf[:,4]-A["kGamma"]*(g-xf[:,1]))+R["g"]*np.cos(g)
    ah=v*np.cos(g)*(xf[:,5]-A["kChi"]*ec)
    command=Rotation.from_euler("ZYX",np.column_stack((chi,g,np.arctan2(ah,ag))))*Rotation.from_euler("y",ac[:,None])
    filt=Rotation.from_quat(x[:,[26,27,28,25]])
    return (command.inv()*filt).magnitude()+(filt.inv()*body).magnitude()

def read_data(group,tag,paths):
    j=load(paths[0])
    if group=="A1":
        s=j["S"];rows=load(paths[1]);Y=s["Y"]
        return s["C"],s["t"],s["x"],ascols(Y,["V","gamma","chi"]),rows["alphaRef"],Y["thetaF"]+Y["thetaR"],Y["kappa"],Y["Tcmd"],rows["gapFeedback"],"completed"
    if group=="WP4":
        d=named(j,"D","names")
        return j["C"],j["t"],j["x"],ascols(d,["actual_V","actual_gamma","actual_chi"]),d["ref_alpha"],d["thetaF"]+d["thetaR"],d["kappa"],d["Traw"],ascols(d,["fb_lower","fb_upper"]),"completed"
    if group=="A3":
        C=load(paths[2])["S"]["C"];t=j["t"];x=j["x"];W=C["MC"]["wind"]
        y,body=air(x,W["mean"]+W["amp"]*np.sin(t[:,None]*W["omega"]+W["phase"]))
        d=named(load(paths[1]),"data","names")
        return C,t,x,y,d["ref_alpha_rad"],attitude(C,x,y,body,d["fb_alpha_rad"]),d["kappa"],d["T_raw_N"],ascols(d,["fb_raw_lower_N","fb_raw_upper_N"]),j["status"]
    C=j["C"];x=j["Xs"];y,body=air(x,0);D=j["D"]
    assert np.max(abs(y[:,0]-D[:,3]))<1e-9
    return C,j["Ts"],x,y,np.deg2rad(D[:,6]),attitude(C,x,y,body,np.deg2rad(D[:,7])),x[:,3]/y[:,0],D[:,10],D[:,14:16],j["Summary"]["status"]

def band(v):
    # Adjacent binary64 values bracket stored-coordinate roundoff. All subsequent
    # algebra is outward fixed-point. This convention does not bound ODE error.
    l=math.nextafter(float(v),-math.inf); h=math.nextafter(float(v),math.inf)
    ln,ld=l.as_integer_ratio();hn,hd=h.as_integer_ratio()
    return I(ln*S//ld,-((-hn*S)//hd))
def status(b,limit):
    return 1 if b.h<=limit.l else (-1 if b.l>limit.h else 0)
def absi(z):return I(0 if z.l<=0<=z.h else min(abs(z.l),abs(z.h)),max(abs(z.l),abs(z.h)))
def writetable(p,rows):
    with p.open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def audit(job):
    group,tag,paths=job; C,t,x,y,refA,theta,kap,raw,gaps,end=read_data(group,tag,paths)
    sol=json.loads((H/"search.json").read_bytes());eps=F(sol["certified_halfwidth"])
    ck=Kernel(eps,sol["domain"]=="subset");match=kernel(C)==MODEL["binary_model"];assert match
    assert C["P7"]["Tmin"]==0 and C["P7"]["Tmax"]==60
    n=len(t);assert np.isfinite(x).all() and np.isfinite(y).all() and np.allclose(t[0],0)
    evalr=wp.Evaluator(C);wp.I.dps=55
    # Inflate by one binary unit so that the older evaluator's binary scale covers
    # the prescribed decimal angle widths too.
    test_eps=np.nextafter(float(eps),np.inf)
    assert F.from_float(float(test_eps))/F.from_float(38.1)>=eps/F("38.1")
    domain=np.zeros(n,dtype=int);box=domain.copy();eligible=np.zeros(n,bool)
    geom=eligible.copy();joint=eligible.copy();JB=eligible.copy()
    old160=eligible.copy();surplus=eligible.copy();baseline_radius=pt(".160");surplus_radius=pt(1)
    attitude_ok=theta<=math.radians(.7);obs=(kap>0)&(raw>0)&(raw<60)
    R=C["R"];xf_all=x[:,19:25];vf,gf,cf,nv,ng,nc=xf_all.T
    mass=R["P"]["m"];grav=R["g"];k=.5*R["rho"]*R["P"]["S"]
    cd=R["P"]["theta"][[0,1,4,7,9]];cl=R["P"]["theta"][[15,16,20,21,22]]
    Xf=mass*(nv+grav*np.sin(gf)); Nf=mass*np.hypot(vf*ng+grav*np.cos(gf),vf*np.cos(gf)*nc)
    def fmap(a):return k*vf**2*(np.polynomial.polynomial.polyval(a,cl)+np.polynomial.polynomial.polyval(a,cd)*np.tan(a))+Xf*np.tan(a)
    refg=np.column_stack((Nf-fmap(R["alphaLo"]),fmap(R["alphaHi"])-Nf))
    Tref=(k*vf**2*np.polynomial.polynomial.polyval(refA,cd)+Xf)/np.cos(refA)
    local_loss=refg-gaps
    use=next(r for r in sol["rows"] if r["repeat"])
    upper=np.array([use["max_lower_pN"],use["max_upper_pN"]])/1e12
    radius=np.zeros(n);ratios=np.zeros(n);Tlo=np.full(n,np.nan);Thi=Tlo.copy()
    e_limit=pt(eps);scale=pt(F.from_float(38.1)).hull(pt("38.1"))
    pi=pt(F("3.14159265358979323846264338327950288419716939937510582097494")).hull(
        pt(F("3.14159265358979323846264338327950288419716939937510582097495")))
    for i,xf in enumerate(x[:,19:25]):
        z=list(map(band,xf));ex=[band(v)-z[j] for j,v in enumerate(y[i])]
        # Select a single local wrap branch and verify the entire interval lies in it.
        turns=round((y[i,2]-xf[2])/(2*math.pi))
        ex[2]=ex[2]-pt(2*turns)*pi
        assert absi(ex[2]).h<pi.l
        dev=[absi(ex[0]),scale*absi(ex[1]),scale*absi(ex[2])]
        width=I(max(e.l for e in dev),max(e.h for e in dev))
        z5=[z[0]-ck.v0,z[1],*z[3:6]]
        dom=ck.domain(z5);domain[i]=status(dom,ck.r);box[i]=status(width,e_limit)
        old160[i]=domain[i]==1 and status(width,baseline_radius)==1
        surplus[i]=domain[i]==1 and status(width,surplus_radius)==1
        if surplus[i]:assert gaps[i,0]>=134.674828-1e-8 and gaps[i,1]>=18.020243-1e-8
        ratios[i]=(dom.l+dom.h)/(ck.r.l+ck.r.h);radius[i]=(width.l+width.h)/(2*S)
        eligible[i]=domain[i]==1 and box[i]==1
        if eligible[i]:
            assert min(gaps[i])>0
            assert np.all(local_loss[i]<=upper+1e-8)
            evalr.prepare(xf,refA[i]);bound=evalr.bounds(test_eps,thrust=True)
            geom[i]=bound["G"];JB[i]=bound["J"]
            joint[i]=JB[i] and attitude_ok[i] and obs[i]
            if JB[i] and attitude_ok[i]:assert obs[i]
            if bound["thrust_lower"] is not None:Tlo[i]=bound["thrust_lower"];Thi[i]=bound["thrust_upper"]
    events={}
    if group=="A1":
        with (A2/"event_comparison.csv").open(encoding="utf-8-sig",newline="") as f:
            events={int(r["index"])-1:r for r in csv.DictReader(f) if int(r["id"])==int(tag[3:]) and int(r["confirmed"])}
    evrows=[]
    for i,old in events.items():
        assert i==round(float(old["time_s"])*20) and abs(t[i]-float(old["time_s"]))<=2*math.ulp(float(t[i]))
        evrows.append(dict(id=int(tag[3:]),time_s=t[i],index=i+1,domain=domain[i],box=box[i],width=radius[i],
            covered=int(eligible[i]),G=int(geom[i]),J=int(joint[i]),attitude_sum_deg=math.degrees(theta[i]),
            kappa=kap[i],Traw_N=raw[i],Tact_N=x[i,35],B1_016=int(old160[i]),KB_extra_1=int(surplus[i]),
            old_local_E=int(old["E_included"]),old_local_J=int(old["J_included"])))
    p=H/"coverage"/f"{group}_{tag}.csv.gz"
    with p.open("wb") as file,gzip.GzipFile(filename="",fileobj=file,mode="wb",mtime=0,compresslevel=3) as gz,io.TextIOWrapper(gz,encoding="ascii",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(["t","domain","error_box","VB_ratio_mid","error_radius_mid","E","G","J_air_thrust_box","attitude_ok","J","theta_sum_deg","kappa","Traw_N","Tact_N","fb_lower_N","fb_upper_N","S4_Tlo_N","S4_Thi_N","Tref_N","lower_loss_N","upper_loss_N","B1_016","KB_extra_1"])
        for i in range(n):w.writerow([t[i],domain[i],box[i],ratios[i],radius[i],int(eligible[i]),int(geom[i]),int(JB[i]),
             int(attitude_ok[i]),int(joint[i]),math.degrees(theta[i]),kap[i],raw[i],x[i,35],*gaps[i],Tlo[i],Thi[i],Tref[i],*local_loss[i],int(old160[i]),int(surplus[i])])
    result=dict(group=group,trajectory=tag,n=n,status=str(end),end_s=t[-1],kernel_match=int(match),
        inside=int(sum(domain==1)),boundary=int(sum(domain==0)),outside=int(sum(domain==-1)),
        error_inside=int(sum(box==1)),error_uncertain=int(sum(box==0)),covered=int(sum(eligible)),G=int(sum(geom)),
        J_air_box=int(sum(JB)),J=int(sum(joint)),covered_attitude_ok=int(sum(eligible&attitude_ok)),
        raw_strict_points=int(sum(obs)),min_kappa=min(kap),min_Traw_N=min(raw),max_Traw_N=max(raw),
        min_Tact_N=min(x[:,35]),max_Tact_N=max(x[:,35]),min_Tref_N=min(Tref),max_Tref_N=max(Tref),
        B1_016=int(sum(old160)),KB_extra_1=int(sum(surplus)),events=len(events),covered_events=sum(r["covered"] for r in evrows))
    return result,evrows

def run():
    (H/"coverage").mkdir(exist_ok=True);jobs=records()
    wp_counts={g:sum(k[0]==g for k in jobs) for g in ["WP4","A1","stress","A3"]}
    assert wp_counts=={"WP4":9,"A1":17,"stress":9,"A3":18},wp_counts
    with ProcessPoolExecutor(12) as executor:vals=list(executor.map(audit,jobs))
    rows=[r for r,ev in vals];events=[e for r,ev in vals for e in ev]
    assert len(events)==21 and all(r["inside"]+r["outside"]+r["boundary"]==r["n"] for r in rows)
    writetable(H/"coverage_by_run.csv",rows);writetable(H/"event_coverage.csv",events)
    print("Classified",len(rows),"records;",len(events),"archived events")
if __name__=="__main__":run()
