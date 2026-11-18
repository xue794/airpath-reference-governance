"""Independent high precision corner checks, identity and manifest checks."""
from pathlib import Path
import sys,json,csv,hashlib,itertools,math
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import mpmath as mp
import numpy as np
from scipy.io import loadmat
mp.mp.dps=70
def m(x):return mp.mpf(float(x))
def poly(c,a):return mp.fsum(v*a**i for i,v in enumerate(c))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run():
    cat=json.loads((ROOT/'inputs/trajectories.json').read_text(encoding='utf-8'))
    manifest=json.loads((ROOT/'inputs/source_manifest.json').read_text(encoding='utf-8'))
    for r in cat:assert sha(r['path'])==r['sha256']
    for r in manifest:assert sha(r['source'])==r['sha256']==sha(r['snapshot'])
    result=dict(input_trajectories_unchanged=len(cat),source_files_unchanged=len(manifest),precision_digits=70,
        sample_rule='12 equally spaced saved indices per trajectory, without selection by result',
        numerical_tolerance=1e-8,checks=[],interval_enclosure_violations=0)
    for rec in cat:
        key=rec['method']+'_'+rec['case'];L=loadmat(ROOT/'results'/f'{key}_export.mat',simplify_cells=True)
        C=L['C'];R=C['R'];A=C['A'];D=L['D']
        with (ROOT/'results'/f'{key}_transfer.csv').open(encoding='utf-8') as f:rr=list(csv.DictReader(f))
        with (ROOT/'results'/f'{key}_bounds.csv').open(encoding='utf-8') as f:bb=list(csv.DictReader(f))
        CL=[m(x) for x in np.array(R['P']['theta'])[[15,16,20,21,22]]];CD=[m(x) for x in np.array(R['P']['theta'])[[0,1,4,7,9]]]
        K=m(.5)*m(R['rho'])*m(R['P']['S']);mass=m(R['P']['m']);g=m(R['g']);scale=m(38.1)
        alo=m(R['alphaLo']);ahi=m(R['alphaHi']);tm=mp.tan(m(R['muMax']));theta=m(.70)*mp.pi/180
        selected=np.unique(np.linspace(0,len(D)-1,12,dtype=int));count=0;ma_count=0;min_distance=mp.inf
        for k in selected:
            xf=[m(v) for v in D[k,1:7]]
            for label in ['E','G','J']:
                eps=m(float(rr[k]['radius_'+label+'_mps']));b=bb[k]
                if eps<=0:continue
                for signs in list(itertools.product([-1,1],repeat=3))+[(0,0,0)]:
                    eV=eps*signs[0];eg=eps*signs[1]/scale;ec=eps*signs[2]/scale
                    V=xf[0]+eV;ga=xf[1]+eg;X=mass*(xf[3]-m(A['kV'])*eV+g*mp.sin(ga))
                    ag=V*(xf[4]-m(A['kGamma'])*eg)+g*mp.cos(ga)
                    ac=V*mp.cos(ga)*(xf[5]-m(A['kChi'])*ec);N=mass*mp.sqrt(ag*ag+ac*ac)
                    def M(a):return K*V*V*(poly(CL,a)+poly(CD,a)*mp.tan(a))+X*mp.tan(a)
                    bm=N-M(alo);bp=M(ahi)-N
                    values=dict(endpoint_lower=bm,endpoint_upper=bp,normal=ag-m(R['aGammaMin']),bank_plus=tm*ag-ac,bank_minus=tm*ag+ac)
                    for name,value in values.items():
                        assert value>=m(float(b[label+'_'+name]))-m(1e-8),(key,k,label,name)
                    for a in [alo,(alo+ahi)/2,ahi]:
                        ma=K*V*V*(poly([i*CL[i] for i in range(1,5)],a)+poly([i*CD[i] for i in range(1,5)],a)*mp.tan(a))+(K*V*V*poly(CD,a)+X)/(mp.cos(a)**2)
                        assert ma>=m(float(b[label+'_Ma']))-m(1e-8)
                        ma_count+=1
                    if label=='J':
                        lo,hi=alo,ahi
                        for _ in range(100):
                            mid=(lo+hi)/2
                            if M(mid)<N:lo=mid
                            else:hi=mid
                        alpha=(lo+hi)/2
                        assert alpha>=m(float(b['J_alpha_lower']))-m(1e-10) and alpha<=m(float(b['J_alpha_upper']))+m(1e-10)
                        num=K*V*V*poly(CD,alpha)+X
                        for kap in [mp.cos(alpha+theta),mp.cos(max(mp.mpf(0),alpha-theta))]:
                            thrust=num/kap
                            assert thrust>=m(float(b['J_thrust_lower']))-m(1e-8) and thrust<=m(float(b['J_thrust_upper']))+m(1e-8)
                    count+=1
        # A separate full-grid pass checks every stored certificate flag/value.
        for row,b in zip(rr,bb):
            for label in ['E','G','J']:
                cert=row['certified_'+label]=='True';r=float(row['radius_'+label+'_mps'])
                assert cert==(r>0)
                if not cert:continue
                assert min(float(b[label+'_endpoint_lower']),float(b[label+'_endpoint_upper']),float(b[label+'_Ma']))>0
                if label in ['G','J']:assert min(float(b[label+'_normal']),float(b[label+'_bank_plus']),float(b[label+'_bank_minus']))>0
                if label=='J':assert float(b['J_thrust_lower'])>0 and float(b['J_thrust_upper'])<60 and float(b['J_kappa_lower'])>0
        result['checks'].append(dict(key=key,selected_indices=selected.tolist(),evaluations=count,derivative_checks=ma_count,all_saved_certificates_checked=len(rr)))
        print('INDEPENDENT_CHECK_PASS',key,count,flush=True)
    result['all_passed']=True
    (ROOT/'checks/independent_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('WP3_INDEPENDENT_VERIFICATION_COMPLETE',flush=True)
if __name__=='__main__':run()
