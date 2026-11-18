"""WP3 conditional instantaneous transfer audit, with outward interval bounds.

No trajectory integration. Binary64 nominal inputs; high-precision mpmath.iv
arithmetic. Not formal verification of the MATLAB simulator or physical plant.
"""
from pathlib import Path
import sys,json,csv,math,time,hashlib,itertools
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import mpmath as mp
import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
mp.iv.dps=40
I=mp.iv
PROTO=json.loads((ROOT/'inputs/protocol.json').read_text(encoding='utf-8'))
SCALE=PROTO['error_box']['speed_scale_mps'];CAP=PROTO['interval']['search_cap_mps'];NBI=PROTO['interval']['bisections']
THETA=PROTO['attitude']['primary_deg']
def point(x):return I.mpf(float(x))
def span(a,b):return I.mpf([float(a),float(b)])
def low(x):return np.nextafter(float(x.a),-np.inf)
def high(x):return np.nextafter(float(x.b),np.inf)
def poly(c,a):
    v=c[-1]
    for x in c[-2::-1]:v=v*a+x
    return v
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
class Evaluator:
    def __init__(self,C):
        R=C['R'];A=C['A'];self.R=R
        self.CL=np.array(R['P']['theta'])[[15,16,20,21,22]];self.CD=np.array(R['P']['theta'])[[0,1,4,7,9]]
        self.cl=[point(x) for x in self.CL];self.cd=[point(x) for x in self.CD]
        self.alo=float(R['alphaLo']);self.ahi=float(R['alphaHi']);self.aa=span(self.alo,self.ahi)
        self.K=point(.5)*point(R['rho'])*point(R['P']['S'])
        self.m=point(R['P']['m']);self.g=point(R['g']);self.kv=point(A['kV']);self.kg=point(A['kGamma']);self.kc=point(A['kChi'])
        self.tm=I.tan(point(R['muMax']));self.agmin=point(R['aGammaMin'])
        self.tanlo=I.tan(point(self.alo));self.tanhi=I.tan(point(self.ahi))
        self.coefflo=poly(self.cl,point(self.alo))+poly(self.cd,point(self.alo))*self.tanlo
        self.coeffhi=poly(self.cl,point(self.ahi))+poly(self.cd,point(self.ahi))*self.tanhi
        self.derivCL=poly([point(j)*self.cl[j] for j in range(1,5)],self.aa)
        self.derivCD=poly([point(j)*self.cd[j] for j in range(1,5)],self.aa)
        self.cdall=poly(self.cd,self.aa);self.tanall=I.tan(self.aa);self.cosall2=I.cos(self.aa)**2
    def numeric(self,xf,e):
        R=self.R;V=xf[0]+e[0];ga=xf[1]+e[1]
        m=float(R['P']['m']);g=float(R['g'])
        X=m*(xf[3]-float(self.kv.mid)*e[0]+g*np.sin(ga))
        ag=V*(xf[4]-float(self.kg.mid)*e[1])+g*np.cos(ga)
        ac=V*np.cos(ga)*(xf[5]-float(self.kc.mid)*e[2]);N=m*np.hypot(ag,ac)
        K=.5*float(R['rho'])*float(R['P']['S']);q=K*V*V
        def M(a):return q*(np.polynomial.polynomial.polyval(a,self.CL)+np.polynomial.polynomial.polyval(a,self.CD)*np.tan(a))+X*np.tan(a)
        return np.array([N-M(self.alo),M(self.ahi)-N]),X,ag,ac,N
    def prepare(self,xf,af):
        self.xf=np.asarray(xf);self.vf,self.gf,self.cf,self.nv,self.ng,self.nc=map(point,xf)
        self.af=point(af);self.af_float=float(af)
        self.cfcoef=poly(self.cl,self.af)+poly(self.cd,self.af)*I.tan(self.af)
        self.tanf=I.tan(self.af)
    def bounds(self,eps,thrust=False,theta_deg=THETA):
        e=span(-eps,eps);eg=e/point(SCALE);ec=eg
        V=self.vf+e;ga=self.gf+eg;q=self.K*V**2
        X=self.m*(self.nv-self.kv*e+self.g*I.sin(ga))
        ag=V*(self.ng-self.kg*eg)+self.g*I.cos(ga)
        ac=V*I.cos(ga)*(self.nc-self.kc*ec)
        N=self.m*I.sqrt(ag**2+ac**2)
        bm=N-q*self.coefflo-X*self.tanlo;bp=q*self.coeffhi+X*self.tanhi-N
        mua=q*self.derivCL+q*self.derivCD*self.tanall+(q*self.cdall+X)/self.cosall2
        b=dict(endpoint_lower=low(bm),endpoint_upper=low(bp),Ma=low(mua),
            normal=low(ag-self.agmin),bank_plus=low(self.tm*ag-ac),bank_minus=low(self.tm*ag+ac),
            V_lower=low(V),gamma_abs_upper=max(abs(low(ga)),abs(high(ga))),
            thrust_lower=None,thrust_upper=None,kappa_lower=None,alpha_lower=None,alpha_upper=None)
        regular=low(V)>0 and b['gamma_abs_upper']<math.pi/2 and b['Ma']>0
        b['E']=regular and min(b['endpoint_lower'],b['endpoint_upper'])>0
        b['G']=b['E'] and min(b['normal'],b['bank_plus'],b['bank_minus'])>0
        b['J']=False
        if not thrust or not b['G']:return b
        residual=q*self.cfcoef+X*self.tanf-N
        mag=max(abs(low(residual)),abs(high(residual)))
        da=high(point(mag)/point(b['Ma']))
        al=max(self.alo,np.nextafter(self.af_float-da,-np.inf));ah=min(self.ahi,np.nextafter(self.af_float+da,np.inf))
        ar=span(al,ah);theta=point(theta_deg)*I.pi/point(180)
        kl=I.cos(point(ah)+theta)
        ul=max(0,np.nextafter(al-high(theta),-np.inf));ku=I.cos(point(ul))
        b.update(alpha_lower=al,alpha_upper=ah,kappa_lower=low(kl))
        if low(kl)<=0 or ah+high(theta)>=math.pi/2:return b
        T=(q*poly(self.cd,ar)+X)/span(low(kl),high(ku))
        b.update(thrust_lower=low(T),thrust_upper=high(T))
        b['J']=low(T)>float(self.R['Tmin']) and high(T)<float(self.R['Tmax'])
        return b
    def radius(self,kind,upper):
        joint=kind=='J';base=self.bounds(0,thrust=joint)
        if not base[kind]:return 0.,False,False,base,'origin_not_certified'
        top=self.bounds(upper,thrust=joint)
        if top[kind]:return upper,upper>0,upper==CAP,top,'passed_upper'
        lo=0.;hi=upper
        for _ in range(NBI):
            mid=(lo+hi)/2;b=self.bounds(mid,thrust=joint)
            if b[kind]:lo=mid
            else:hi=mid
        final=self.bounds(lo,thrust=joint)
        assert final[kind]
        return lo,lo>0,False,final,'positive' if lo>0 else 'below_search_resolution'

def analyze(key):
    start=time.time();f=ROOT/'results'/f'{key}_export.mat';L=loadmat(f,simplify_cells=True)
    C=L['C'];E=Evaluator(C);D=L['D'];names=[str(x) for x in L['names']];idx={v:i for i,v in enumerate(names)}
    xs=L['x'];rows=[];boundrows=[];ident_max=np.zeros(9)
    mode=key.split('_')[-1]
    for k,row in enumerate(D):
        d=dict(zip(names,row));xf=row[1:7];t=d['time_s'];x=xs[k]
        # Independent wind-aware air-state reconstruction from saved plant state.
        rot=Rotation.from_quat(np.r_[x[7:10],x[6]])
        wind=np.zeros(3) if mode=='hold' else np.array([6+2*np.sin(.25*t),-3+1.5*np.sin(.40*t+.5),.8*np.sin(.2*t)])
        vaI=rot.apply(x[3:6])-wind;vaB=rot.inv().apply(vaI);V=np.linalg.norm(vaI)
        gamma=np.arctan2(-vaI[2],np.hypot(vaI[0],vaI[1]));chi=np.arctan2(vaI[1],vaI[0]);kap=vaB[0]/V
        e=np.array([V-xf[0],gamma-xf[1],np.arctan2(np.sin(chi-xf[2]),np.cos(chi-xf[2]))])
        hf,Xf,agf,acf,Nf=E.numeric(xf,np.zeros(3));hc,Xc,agc,acc,Nc=E.numeric(xf,e)
        K=.5*float(E.R['rho'])*float(E.R['P']['S']);angles=np.array([E.alo,E.ahi])
        ce=np.polynomial.polynomial.polyval(angles,E.CL)+np.polynomial.polynomial.polyval(angles,E.CD)*np.tan(angles)
        dm=K*(2*xf[0]*e[0]+e[0]**2)*ce+(Xc-Xf)*np.tan(angles);dn=Nc-Nf
        identity=hf+np.array([dn-dm[0],dm[1]-dn])
        cd=np.polynomial.polynomial.polyval(d['alphaC'],E.CD)
        Traw=(K*V*V*cd+Xc)/kap
        errs=[np.max(abs(vaI-rot.apply(x[3:6])+wind)),np.max(abs(wind-row[40:43])),abs(V-d['actual_V']),
            np.max(abs(e-row[7:10])),np.max(abs(hf-row[17:19])),np.max(abs(hc-row[26:28])),np.max(abs(identity-hc)),abs(Traw-d['Traw']),abs(kap-d['kappa'])]
        ident_max=np.maximum(ident_max,errs)
        assert max(errs)<1e-7,(key,k,errs)
        E.prepare(xf,d['ref_alpha'])
        rE,ceE,capE,bE,sE=E.radius('E',CAP)
        rG,ceG,capG,bG,sG=E.radius('G',rE)
        rJ,ceJ,capJ,bJ,sJ=E.radius('J',rG)
        assert rJ<=rG<=rE
        eb=max(abs(e[0]),SCALE*abs(e[1]),SCALE*abs(e[2]));tc=d['thetaC']*180/math.pi;ts=(d['thetaF']+d['thetaR'])*180/math.pi
        inE=ceE and eb<=rE;inG=ceG and eb<=rG;airJ=ceJ and eb<=rJ
        sumJ=airJ and ts<=THETA;directJ=airJ and tc<=THETA
        observed=E.bounds(np.nextafter(eb,np.inf),thrust=True,theta_deg=np.nextafter(tc,np.inf))
        if inE:assert min(hc)>-1e-8
        if inG:assert min(agc-float(E.R['aGammaMin']),np.tan(float(E.R['muMax']))*agc-abs(acc))>-1e-8
        if directJ:assert -1e-8<d['Traw']<60+1e-8
        rr=dict(time_s=t,error_box_mps=eb,theta_sum_deg=ts,theta_direct_deg=tc,
            radius_E_mps=rE,radius_G_mps=rG,radius_J_mps=rJ,certified_E=ceE,certified_G=ceG,certified_J=ceJ,
            capped_E=capE,capped_G=capG,capped_J=capJ,inclusion_E=inE,inclusion_G=inG,inclusion_J_air=airJ,
            inclusion_J_S4=sumJ,inclusion_J_direct=directJ,attitude_S4_ok=ts<=THETA,attitude_direct_ok=tc<=THETA,
            observed_box_joint_pass=bool(observed['J']),status_E=sE,status_G=sG,status_J=sJ,
            lower_loss_N=hf[0]-hc[0],upper_loss_N=hf[1]-hc[1],deltaN_N=dn,deltaMlo_N=dm[0],deltaMhi_N=dm[1],
            ref_lower_N=hf[0],ref_upper_N=hf[1],fb_lower_N=hc[0],fb_upper_N=hc[1],Traw_N=d['Traw'],Tstar_N=d['Tstar'])
        rows.append(rr)
        bb={'time_s':t}
        for label,bnd in [('E',bE),('G',bG),('J',bJ)]:
            for n,v in bnd.items():bb[label+'_'+n]=v
        boundrows.append(bb)
        if k%400==0:print(f'{key} {k}/{len(D)} elapsed {time.time()-start:.1f}s',flush=True)
    for suffix,data in [('transfer',rows),('bounds',boundrows)]:
        with (ROOT/'results'/f'{key}_{suffix}.csv').open('w',newline='',encoding='utf-8') as out:
            w=csv.DictWriter(out,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
    def a(n):return np.array([r[n] for r in rows])
    result=dict(key=key,case=mode,method=key.rsplit('_',1)[0],samples=len(rows),end_s=float(D[-1,0]),elapsed_s=time.time()-start,
        reconstruction_max_errors=ident_max.tolist(),radius_stats={},counts={},max_error_box_mps=float(a('error_box_mps').max()),
        max_theta_sum_deg=float(a('theta_sum_deg').max()),max_theta_direct_deg=float(a('theta_direct_deg').max()))
    for label in ['E','G','J']:
        v=a('radius_'+label+'_mps');result['radius_stats'][label]=dict(min=float(v.min()),median=float(np.median(v)),max=float(v.max()))
    for n in rows[0]:
        if n.startswith(('certified','capped','inclusion','attitude','observed')):result['counts'][n]=int(a(n).sum())
    result['failure_attribution_J_S4']={
        'uncertified':int((~a('certified_J')).sum()),
        'air_only':int((a('certified_J')&~a('inclusion_J_air')&a('attitude_S4_ok')).sum()),
        'attitude_only':int((a('inclusion_J_air')&~a('attitude_S4_ok')).sum()),
        'both':int((a('certified_J')&~a('inclusion_J_air')&~a('attitude_S4_ok')).sum()),
        'covered':int(a('inclusion_J_S4').sum())}
    assert sum(result['failure_attribution_J_S4'].values())==len(rows)
    dump(ROOT/'results'/f'{key}_summary.json',result)
    print(f'COMPLETE {key}: {result["counts"]}',flush=True)
    return result
if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    cat=json.loads((ROOT/'inputs/trajectories.json').read_text(encoding='utf-8'))
    keys=[r['method']+'_'+r['case'] for r in cat]
    if '--key' in sys.argv:keys=[sys.argv[sys.argv.index('--key')+1]]
    with ProcessPoolExecutor(max_workers=3) as pool:out=list(pool.map(analyze,keys))
    dump(ROOT/'results/transfer_summary.json',out)
    print('WP3_INTERVAL_COMPLETE',flush=True)
