"""Targeted independent check of the newly observed negative reference thrust."""
from pathlib import Path
import sys,json,csv
import numpy as np
from scipy.io import loadmat
H=Path(__file__).resolve().parent;sys.path.insert(0,str(H/'analysis/vendor'))
import mpmath as mp
mp.mp.dps=70
L=loadmat(H/'analysis/results/Transparent_noisy_export.mat',simplify_cells=True)
C=L['C'];R=C['R'];D=L['D'];ix={str(n):i for i,n in enumerate(L['names'])}
with (H/'analysis/results/Transparent_noisy_transfer.csv').open() as f:tr=list(csv.DictReader(f))
neg=D[:,ix['ref_T']]<-1e-8;zero=np.array([float(r['radius_J_mps'])==0 for r in tr])
assert np.array_equal(neg,zero)
ids=np.flatnonzero(neg);selected=np.unique(np.clip(np.r_[ids,ids[0]-1,ids[-1]+1],0,len(D)-1))
def m(x):return mp.mpf(float(x))
def poly(cc,a):return mp.fsum(c*a**i for i,c in enumerate(cc))
CL=[m(x) for x in np.array(R['P']['theta'])[[15,16,20,21,22]]];CD=[m(x) for x in np.array(R['P']['theta'])[[0,1,4,7,9]]]
K=m(.5)*m(R['rho'])*m(R['P']['S']);mass=m(R['P']['m']);g=m(R['g']);rows=[]
for k in selected:
    V,ga,chi,nv,ng,nc=map(m,D[k,1:7]);X=mass*(nv+g*mp.sin(ga));ag=V*ng+g*mp.cos(ga);ac=V*mp.cos(ga)*nc;N=mass*mp.sqrt(ag*ag+ac*ac);q=K*V*V
    def M(a):return q*(poly(CL,a)+poly(CD,a)*mp.tan(a))+X*mp.tan(a)
    lo=m(R['alphaLo']);hi=m(R['alphaHi']);assert M(lo)<N<M(hi)
    for _ in range(150):
        mid=(lo+hi)/2
        if M(mid)<N:lo=mid
        else:hi=mid
    a=(lo+hi)/2;T=(q*poly(CD,a)+X)/mp.cos(a);diff=abs(float(T)-D[k,ix['ref_T']]);assert diff<1e-10
    assert (T<0)==bool(neg[k])
    rows.append(dict(index=int(k),time_s=float(D[k,0]),reference_T_N=float(T),stored_reference_T_N=float(D[k,ix['ref_T']]),absolute_difference_N=diff,
        reference_alpha_deg=float(a*180/mp.pi),feedback_raw_T_N=float(D[k,ix['Traw']]),joint_radius_mps=float(tr[k]['radius_J_mps'])))
out=dict(precision_digits=70,selection='all 17 negative saved points plus the adjacent saved points; follow-up diagnostic of a new finding',
    negative_samples=int(neg.sum()),first_negative_sample_s=float(D[ids[0],0]),last_negative_sample_s=float(D[ids[-1],0]),
    zero_radius_indices_equal_negative_reference_indices=True,max_difference_N=max(r['absolute_difference_N'] for r in rows),
    all_checked_roots_inside_original_interval=True,rows=rows,continuous_event_times_claimed=False)
(H/'checks/negative_reference_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='rows'}))
