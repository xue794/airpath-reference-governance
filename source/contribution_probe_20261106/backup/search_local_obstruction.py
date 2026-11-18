"""Optimizer proposes points; original outward interval kernel verifies signs.
No optimizer result by itself is labeled a certificate.
"""
from pathlib import Path
import json,sys
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize
from scipy.io import loadmat
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'original_certificate'))
import interval_certificate as ic
R=loadmat(HERE/'V31_backup_saved_audit.mat',simplify_cells=True)['R']
P=R['PB'];L=np.linalg.cholesky(np.linalg.inv(P));m=R['P']['m'];g=R['g']
C=.5*R['rho']*R['P']['S'];a=R['alphaHi'];ta=np.tan(a)
CD=.019+a*(-.078+a*(3.420+a*(1.819-6.563*a)))
CL=.016+a*(5.343+a*(.506+a*(-36.30+46.13*a)))

def values(z,center=31):
    dv,ga,nv,ng,nc=z;V=center+dv;cg=np.cos(ga);sg=np.sin(ga)
    uv=-.6*dv-1.4*nv;ug=-1.2*ga-1.8*ng;uc=-1.5*nc
    ag=V*ng+g*cg;ac=V*cg*nc;agd=nv*ng+V*ug-g*sg*ng
    acd=(nv*cg-V*sg*ng)*nc+V*cg*uc
    X=m*(nv+g*sg);Xd=m*(uv+g*cg*ng);rr=np.hypot(ag,ac)
    N=m*rr;Nd=m*(ag*agd+ac*acd)/rr
    D=C*V*V*CD;lift=C*V*V*CL;M=lift+(D+X)*ta
    h=M-N-10;phi=2*(lift+D*ta)/V*nv+ta*Xd-Nd+h
    return h,phi

results=[]
for rho in [.005,.01]:
    for target,label in [(0,'hNplus'),(1,'phiNplus')]:
        proposals=[]
        for seed in [np.zeros(5),np.array([0,0,0,-.9,0]),np.array([0,0,0,.9,0]),np.array([0,-.9,0,0,0])]:
            sol=minimize(lambda y:values(np.sqrt(rho)*L@y)[target],seed,method='SLSQP',
                constraints={'type':'ineq','fun':lambda y:.999**2-y@y},
                options={'ftol':1e-11,'maxiter':250})
            proposals.append(sol)
        sol=min(proposals,key=lambda v:v.fun)
        z=np.sqrt(rho)*L@sol.x
        ic.Vb=ic.I('31'); zi=[ic.I(str(Fraction.from_float(float(v)))) for v in z]
        bounds,valid=ic.evaluate(zi);ell=ic.ellipse(zi)
        negative=label in bounds and bounds[label].hi<0 and ell.hi<ic.I(str(rho)).lo
        results.append(dict(center=31,rho=rho,target=label,optimizer_success=bool(sol.success),
            optimizer_message=str(sol.message),proposed_value=float(sol.fun),xi=z.tolist(),
            inside_ellipsoid_interval=ell.pair(),root_and_10N_state_reserve_checked=bool(valid),
            negative_sign_and_interior_strictly_verified=bool(negative),
            target_interval=bounds.get(label).pair(),
            all_checks={k:v.pair() for k,v in bounds.items()}))
        print(rho,label,float(sol.fun),'strict negative',negative,flush=True)
(HERE/'local_obstructions.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
