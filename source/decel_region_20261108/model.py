"""Frozen inputs, numerical scouting, and exact algebraic tube geometry."""
import hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq,minimize
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OLD=ROOT/'decel_screening_20261108'
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
data=loadmat(ROOT/'endpoint_probe_20261106/full6dof/V31_C60_warm_fine/N10.mat',simplify_cells=True)
R=data['C']['R']; A=data['C']['A']; P=R['PB'];m=R['P']['m'];g=R['g'];K=.5*R['rho']*R['P']['S']
cl=R['P']['theta'][[15,16,20,21,22]];cd=R['P']['theta'][[0,1,4,7,9]]
def save(p,value):p.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def lock():
    sys.path.insert(0,str(OLD)); from freeze import lock as oldlock
    oldlock()
    for name,h in json.loads((OLD/'release_manifest.json').read_bytes()).items():
        assert hashlib.sha256((OLD/name).read_bytes()).hexdigest()==h
    plan=dict(duration=22,sequence=[.01,.003],fallback=[26,.01,.003],stages=[[0,3],[3,25],[25,35],[35,55]],window=[14,20],
              maxnodes=1000000,cpu_limit_s=14400,endpoint_N=10,thrust_range_N=[3,57],error_box=[.125,.45,.35],
              previous_manifest=hashlib.sha256((OLD/'release_manifest.json').read_bytes()).hexdigest(),
              method='Seek verified interior counterexamples; then outward branch and bound. Stop at the first fully passing candidate. No flight simulation.')
    p=HERE/'protocol.json'
    if p.exists():assert json.loads(p.read_bytes())==plan
    else:save(p,plan)
def center(t,D):
    a,b,h,q=R['Vb'],31.25,D,np.clip((t-3)/D,0,1)
    if t>=35:a,b,h,q=31.25,R['Vb'],20,np.clip((t-35)/20,0,1)
    return a+(b-a)*q**3*(10-15*q+6*q*q),(b-a)*30*q*q*(1-q)**2/h,(b-a)*60*q*(1-q)*(1-2*q)/h**2
def scout(D,r):
    L=np.linalg.cholesky(P).T; M=np.linalg.inv(L)*np.sqrt(r)*.9999
    def force(v,ga,nv,ng,nc):
        x=m*(nv+g*np.sin(ga));n=m*np.hypot(v*ng+g*np.cos(ga),v*np.cos(ga)*nc)
        def lift(a):return K*v*v*np.polynomial.polynomial.polyval(a,cl)
        def drag(a):return K*v*v*np.polynomial.polynomial.polyval(a,cd)
        def gap(a):return lift(a)+(drag(a)+x)*np.tan(a)-n
        return gap,drag,x
    def values(y):
        v,nv,_=center(y[0],D);dv,ga,dn,ng,nc=M@y[1:]
        gap,drag,x=force(v+dv,ga,nv+dn,ng,nc); alpha=brentq(gap,-.05,.22)
        return (drag(alpha)+x)/np.cos(alpha)-3,gap(R['alphaHi'])-10
    for j,start in [(0,3+D/2),(1,3+D)]:
        res=minimize(lambda y:values(y)[j],np.r_[start,np.zeros(5)],method='SLSQP',bounds=[(3,35)]+[(-1,1)]*5,
                     constraints={'type':'ineq','fun':lambda y:1-y[1:]@y[1:]},options={'ftol':1e-10,'maxiter':100})
        if values(res.x)[j]<-1e-6:
            e=M@res.x[1:]; assert e@P@e<r
            return dict(D=D,r=r,t=float(res.x[0]),e=e.tolist(),item=['hTmin','hNplus'][j],scout=values(res.x)[j])
def geometry(r):
    p=[[F(float(v)) for v in row] for row in P]
    B=[[F(float(v)) for v in row] for row in R['AB']]
    Q=[[-sum(B[k][i]*p[k][j]+p[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
    assert all(p[i][j]==p[j][i] and Q[i][j]==Q[j][i] for i in range(5) for j in range(5))
    def gersh(a):return min(row[i]-sum(abs(v) for j,v in enumerate(row) if j!=i) for i,row in enumerate(a))
    # Directly verifies -A'P-PA >= 0.12 P for stored binary64 coefficients.
    H=[[Q[i][j]-F('0.12')*p[i][j] for j in range(5)] for i in range(5)]
    assert gersh(H)>0 and gersh(p)>0 and r<=R['rhoB'] and np.all(data['Xs'][0,19:25]==[R['Vb'],0,0,0,0,0])
    assert all(p[i][j]==0 for i in range(5) for j in range(5) if i!=j and {i,j} not in [{0,2},{1,3}])
    inv=[None]*5
    for i,j in [(0,2),(1,3)]:
        det=p[i][i]*p[j][j]-p[i][j]**2
        inv[i]=p[j][j]/det;inv[j]=p[i][i]/det
    inv[4]=1/p[4][4]
    half=[np.nextafter(np.sqrt(float(F(r)*v)),np.inf) for v in inv]
    assert all(F(h)**2>=F(r)*v for h,v in zip(half,inv))
    return dict(projection=half,lyapunov_rate=.12,gershgorin_H=str(gersh(H)),terminal_slack=float(R['rhoB']-r),initial_error=0)
if __name__=='__main__':
    lock()
    for r in [.01,.003]:
        w=scout(22,r); print(r,w,geometry(r),flush=True)
        if w:save(OUT/('scout_'+str(r)+'.json'),w)

