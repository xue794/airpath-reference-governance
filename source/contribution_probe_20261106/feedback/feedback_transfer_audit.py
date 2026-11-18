"""Exact margin-loss identities and interval-enclosed local transfer tubes.

All inputs are stored second-round trajectories. No aircraft integration.
mpmath.iv operates on fixed nominal binary floating-point input constants.
This is a numerical interval calculation, not a formally verified software proof.
"""
from pathlib import Path
import sys,json,csv,math,itertools,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import mpmath as mp
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation
mp.iv.dps=40
I=mp.iv
SOURCE=ROOT.parents[1]/'endpoint_probe_20261106'/'full6dof'
SPEED_SCALE=38.1
THETA_DEG=.70

def scalar(v):return I.mpf([float(v),float(v)])
def span(a,b):return I.mpf([float(a),float(b)])
def low(v):return np.nextafter(float(v.a),-np.inf)
def high(v):return np.nextafter(float(v.b),np.inf)
def horner(c,x):
    z=scalar(c[-1])
    for a in c[-2::-1]:z=z*x+scalar(a)
    return z
def poly(c,x):return np.polynomial.polynomial.polyval(x,c)
def prepare(C):
    R=C['R'];th=R['P']['theta'];A=C['A']
    return dict(m=float(R['P']['m']),g=float(R['g']),K=.5*float(R['rho'])*float(R['P']['S']),
        CL=np.array(th)[[15,16,20,21,22]],CD=np.array(th)[[0,1,4,7,9]],
        alo=float(R['alphaLo']),ahi=float(R['alphaHi']),kV=float(A['kV']),kG=float(A['kGamma']),kC=float(A['kChi']))
def forces(V,ga,eg,ec,xf,P):
    X=P['m']*(xf[3]-P['kV']*(V-xf[0])+P['g']*np.sin(ga))
    ag=V*(xf[4]-P['kG']*eg)+P['g']*np.cos(ga)
    ac=V*np.cos(ga)*(xf[5]-P['kC']*ec)
    return X,ag,ac,P['m']*np.hypot(ag,ac)
def mmap(a,V,X,P):return P['K']*V*V*(poly(P['CL'],a)+poly(P['CD'],a)*np.tan(a))+X*np.tan(a)
def root_alpha(V,X,N,P):
    flo=mmap(P['alo'],V,X,P)-N;fhi=mmap(P['ahi'],V,X,P)-N
    if flo>1e-6 or fhi<-1e-6:return float('nan')
    if abs(flo)<1e-9:return P['alo']
    if abs(fhi)<1e-9:return P['ahi']
    return brentq(lambda a:mmap(a,V,X,P)-N,P['alo'],P['ahi'],xtol=1e-14)
def calc(xf,e,P):
    V=xf[0]+e[0];ga=xf[1]+e[1];X,ag,ac,N=forces(V,ga,e[1],e[2],xf,P)
    h=np.array([N-mmap(P['alo'],V,X,P),mmap(P['ahi'],V,X,P)-N])
    return h,X,ag,ac,N
def iv_bounds(xf,P,eps,theta_deg=THETA_DEG,include_thrust=True):
    """|eV|<=eps, |egamma|,|echi|<=eps/38.1, dSO3(Cc,C)<=theta.
    Returns enclosures of feedback raw endpoint margins and raw physical thrust.
    The root is bounded by a uniform derivative lower bound and interval residual.
    """
    eV=span(-eps,eps);eg=eV/scalar(SPEED_SCALE);ec=eg
    V=scalar(xf[0])+eV;gam=scalar(xf[1])+eg
    mass=scalar(P['m']);g=scalar(P['g']);qS=scalar(P['K'])*V**2
    X=mass*(scalar(xf[3])-scalar(P['kV'])*eV+g*I.sin(gam))
    ag=V*(scalar(xf[4])-scalar(P['kG'])*eg)+g*I.cos(gam)
    ac=V*I.cos(gam)*(scalar(xf[5])-scalar(P['kC'])*ec)
    N=mass*I.sqrt(ag**2+ac**2)
    def M(a):return qS*(horner(P['CL'],a)+horner(P['CD'],a)*I.tan(a))+X*I.tan(a)
    hm=N-M(scalar(P['alo']));hp=M(scalar(P['ahi']))-N
    raw=(low(hm),low(hp))
    if not include_thrust:return raw+(float('nan'),float('nan'))
    if min(raw)<=0:return raw+(float('-inf'),float('inf'))
    aa=span(P['alo'],P['ahi'])
    cld=sum(scalar(i)*scalar(P['CL'][i])*aa**(i-1) for i in range(1,5))
    cdd=sum(scalar(i)*scalar(P['CD'][i])*aa**(i-1) for i in range(1,5))
    Ma=qS*cld+qS*cdd*I.tan(aa)+(qS*horner(P['CD'],aa)+X)/(I.cos(aa)**2)
    mu=low(Ma)
    if mu<=0:return raw+(float('-inf'),float('inf'))
    hf,Xf,_,_,Nf=calc(xf,np.zeros(3),P);af=root_alpha(xf[0],Xf,Nf,P)
    if not np.isfinite(af):return raw+(float('-inf'),float('inf'))
    res=M(scalar(af))-N
    # Residual/monotonicity bound; outward rounding is kept in the division.
    mag=max(abs(low(res)),abs(high(res)))
    da=high(scalar(mag)/scalar(mu))
    al=max(P['alo'],np.nextafter(af-da,-np.inf));ah=min(P['ahi'],np.nextafter(af+da,np.inf))
    ar=span(al,ah)
    numerator=qS*horner(P['CD'],ar)+X
    theta=scalar(theta_deg)*I.pi/180
    # Rotation distance bounds angle between body x-axis and actual air direction.
    # |acos(kappa)-alpha_c|<=theta; stronger than a generic cosine Lipschitz bound.
    kl=I.cos(scalar(ah)+theta)
    upper_angle=max(0,np.nextafter(al-high(theta),-np.inf));ku=I.cos(scalar(upper_angle))
    if low(kl)<=0:return raw+(float('-inf'),float('inf'))
    kap=span(low(kl),high(ku));Traw=numerator/kap
    return raw+(low(Traw),high(Traw))
def radius(xf,P,kind,theta_deg=THETA_DEG):
    thrust=kind!='endpoints'
    def feasible(e):
        a,b,c,d=iv_bounds(xf,P,e,theta_deg=theta_deg,include_thrust=thrust)
        if a<=0 or b<=0:return False
        if not thrust:return True
        return c>=(3 if kind=='reserve' else 0) and d<=(57 if kind=='reserve' else 60)
    if not feasible(0):return 0.
    lo=0.;hi=1.
    if feasible(hi):return hi
    for _ in range(20):
        mid=(lo+hi)/2
        if feasible(mid):lo=mid
        else:hi=mid
    return lo
def actual_air(x):
    rot=Rotation.from_quat(np.r_[x[7:10],x[6]])
    body=x[3:6];inertial=rot.apply(body);V=np.linalg.norm(body)
    gam=np.arctan2(-inertial[2],np.hypot(*inertial[:2]));chi=np.arctan2(inertial[1],inertial[0]);kap=body[0]/V
    return V,gam,chi,kap,rot
def Raxis(axis,a):return Rotation.from_rotvec(np.eye(3)[axis]*a)
def analyze_case(args):
    case,method=args;path=SOURCE/case/(method+'.mat');L=loadmat(path,simplify_cells=True);P=prepare(L['C']);ts=L['Ts'];xs=L['Xs'];saved=L['D']
    rows=[];max_id_error=0.;max_saved_error=0.;n=len(ts)
    out=ROOT/f'{case}_{method}_transfer.csv'
    summarize_existing='--summarize-existing' in sys.argv
    if summarize_existing:
        with out.open(newline='',encoding='utf-8') as f:
            rows=[{k:(v=='True' if k.endswith('_inclusion') else float(v)) for k,v in r.items()} for r in csv.DictReader(f)]
        for k,r in enumerate(rows):
            r['endpoint_tube_inclusion']=r['radius_endpoints_mps']>0 and r['error_box_mps']<=r['radius_endpoints_mps']
            r['physical_tube_inclusion']=r['radius_four_physical_mps']>0 and r['error_box_mps']<=r['radius_four_physical_mps'] and r['attitude_error_deg']<=THETA_DEG
            r['reserved_tube_inclusion']=r['radius_four_reserved_mps']>0 and r['error_box_mps']<=r['radius_four_reserved_mps'] and r['attitude_error_deg']<=THETA_DEG
            predicted=np.array([r['reference_lower_N']+r['delta_N_N']-r['delta_M_lower_N'],r['reference_upper_N']+r['delta_M_upper_N']-r['delta_N_N']])
            hc=np.array([r['feedback_lower_N'],r['feedback_upper_N']])
            max_id_error=max(max_id_error,float(np.max(abs(predicted-hc))))
            max_saved_error=max(max_saved_error,float(np.max(abs(np.array([r['reference_lower_N'],r['reference_upper_N'],*hc])-saved[k,12:16]))))
    for k,(t,x) in enumerate(zip(ts,xs)):
        if summarize_existing:break
        xf=x[19:25];Va,ga,chi,kap,rot=actual_air(x)
        e=np.array([Va-xf[0],ga-xf[1],np.arctan2(np.sin(chi-xf[2]),np.cos(chi-xf[2]))])
        hf,Xf,agf,acf,Nf=calc(xf,np.zeros(3),P);hc,Xc,agc,acc,Nc=calc(xf,e,P)
        deltaX=Xc-Xf;deltaN=Nc-Nf
        aa=np.array([P['alo'],P['ahi']]);Ce=poly(P['CL'],aa)+poly(P['CD'],aa)*np.tan(aa)
        deltaM=P['K']*(2*xf[0]*e[0]+e[0]**2)*Ce+deltaX*np.tan(aa)
        predicted=hf+np.array([deltaN-deltaM[0],deltaM[1]-deltaN])
        max_id_error=max(max_id_error,float(np.max(abs(predicted-hc))))
        max_saved_error=max(max_saved_error,float(np.max(abs(np.r_[hf,hc]-saved[k,12:16]))))
        af=root_alpha(xf[0],Xf,Nf,P);ac=root_alpha(Va,Xc,Nc,P)
        Tc=P['K']*Va*Va*poly(P['CD'],ac)+Xc if np.isfinite(ac) else np.nan
        Tstar=Tc/np.cos(ac) if np.isfinite(ac) else np.nan;Traw=Tc/kap
        mu=np.arctan2(acc,agc)
        command=Raxis(2,chi)*Raxis(1,ga)*Raxis(0,mu)*Raxis(1,ac)
        theta=(command.inv()*rot).magnitude()*180/np.pi
        errbox=max(abs(e[0]),SPEED_SCALE*abs(e[1]),SPEED_SCALE*abs(e[2]))
        rE=radius(xf,P,'endpoints');rP=radius(xf,P,'physical');rR=radius(xf,P,'reserve')
        # Radius tests use whole error boxes. Direct reconstructed identities are
        # separately checked against saved MATLAB data; inclusion is checked below.
        row=dict(time_s=float(t),eV_mps=e[0],egamma_deg=e[1]*180/np.pi,echi_deg=e[2]*180/np.pi,error_box_mps=errbox,
            attitude_error_deg=theta,reference_lower_N=hf[0],reference_upper_N=hf[1],feedback_lower_N=hc[0],feedback_upper_N=hc[1],
            loss_lower_N=hf[0]-hc[0],loss_upper_N=hf[1]-hc[1],delta_X_N=deltaX,delta_N_N=deltaN,delta_M_lower_N=deltaM[0],delta_M_upper_N=deltaM[1],
            Tstar_N=Tstar,Traw_N=Traw,radius_endpoints_mps=rE,radius_four_physical_mps=rP,radius_four_reserved_mps=rR,
            endpoint_tube_inclusion=rE>0 and errbox<=rE,physical_tube_inclusion=rP>0 and errbox<=rP and theta<=THETA_DEG,
            reserved_tube_inclusion=rR>0 and errbox<=rR and theta<=THETA_DEG)
        rows.append(row)
        if k%100==0:print(f'{case} {method} {k}/{n}',flush=True)
    with out.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    def vals(key):return np.array([r[key] for r in rows])
    summary=dict(case=case,method=method,samples=n,end_time_s=float(ts[-1]),max_exact_identity_error_N=max_id_error,max_saved_MATLAB_margin_error_N=max_saved_error,
        min_endpoint_radius_mps=float(vals('radius_endpoints_mps').min()),max_endpoint_radius_mps=float(vals('radius_endpoints_mps').max()),
        min_four_physical_radius_mps=float(vals('radius_four_physical_mps').min()),min_four_reserved_radius_mps=float(vals('radius_four_reserved_mps').min()),
        endpoint_tube_inclusion_samples=int(vals('endpoint_tube_inclusion').sum()),physical_tube_inclusion_samples=int(vals('physical_tube_inclusion').sum()),reserved_tube_inclusion_samples=int(vals('reserved_tube_inclusion').sum()),
        max_error_box_mps=float(vals('error_box_mps').max()),max_attitude_error_deg=float(vals('attitude_error_deg').max()),
        min_feedback_upper_N=float(vals('feedback_upper_N').min()),max_upper_margin_loss_N=float(vals('loss_upper_N').max()),
        max_lower_margin_loss_N=float(vals('loss_lower_N').max()),min_Traw_N=float(vals('Traw_N').min()),max_Traw_N=float(vals('Traw_N').max()),
        raw_thrust_reserved_violation_samples=int(((vals('Traw_N')<3-1e-6)|(vals('Traw_N')>57+1e-6)).sum()),
        raw_thrust_hardware_violation_samples=int(((vals('Traw_N')<0-1e-6)|(vals('Traw_N')>60+1e-6)).sum()),
        last_sample=rows[-1],input_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    # Validate that a positive interval tube never claims inclusion at a negative endpoint.
    for r in rows:
        if r['endpoint_tube_inclusion']:assert min(r['feedback_lower_N'],r['feedback_upper_N'])>=-1e-6
        if r['physical_tube_inclusion']:assert -1e-6<=r['Traw_N']<=60+1e-6
    (ROOT/f'{case}_{method}_summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item()),encoding='utf-8')
    return summary

if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    cases=[(c,m) for c in ['V31_C60_warm_fine','V32_C60_warm_fine'] for m in ['T','N0','N10']]
    with ProcessPoolExecutor(max_workers=3) as pool:results=list(pool.map(analyze_case,cases))
    (ROOT/'transfer_summary.json').write_text(json.dumps(results,indent=2,default=lambda x:x.item()),encoding='utf-8')
    print('FEEDBACK_TRANSFER_AUDIT_COMPLETE',flush=True)
