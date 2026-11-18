"""Numerically locate candidates; only independently enclosed in-domain points count."""
from fractions import Fraction as F
from itertools import product
import json,math,os
# The tiny numerical candidate searches use one BLAS thread for reproducible
# directions; the two outward interval checks are independent of BLAS.
os.environ["OPENBLAS_NUM_THREADS"]="1"
import numpy as np
from scipy.optimize import minimize
from boxes import H,MODEL,Kernel,pt,S
from independent import Check,q,fraction
from search import writej

def run():
    data=json.loads((H/"search.json").read_bytes())
    b=MODEL["binary_model"]; num=lambda k:float.fromhex(b[k])
    P=np.array([[float.fromhex(v) for v in row] for row in b["PB"]])
    Q=np.linalg.solve(np.linalg.cholesky(P).T,np.eye(5))*math.sqrt(num("rhoB"))
    cl=list(map(float.fromhex,b["CL"]));cd=list(map(float.fromhex,b["CD"]))
    a=num("alphaHi");ta=math.tan(a)
    mass=num("mass");g=num("g");K=num("density")*num("S")/2
    B=np.polynomial.polynomial.polyval(a,cl)+ta*np.polynomial.polynomial.polyval(a,cd)
    def gaps(z,e):
        dv,gam,nv,ng,nc=z;ev,eg,ec=e;v=num("Vb")+dv+ev;gam+=eg
        X=mass*(nv-num("kV")*ev+g*math.sin(gam))
        N=mass*math.hypot(v*(ng-num("kGamma")*eg)+g*math.cos(gam),v*math.cos(gam)*(nc-num("kChi")*ec))
        return np.array([N-K*v*v*cl[0],K*v*v*B+X*ta-N])
    # Derive reference directions for each signed error corner and both endpoint losses.
    templates=[]
    for sign in product([-1,1],repeat=3):
        e=np.array(sign)*np.array([1,1/38.1,1/38.1])
        for side in [0,1]:
            def objective(y):
                z=Q@y
                return -(gaps(z,[0,0,0])-gaps(z,.2*e))[side]
            res=minimize(objective,np.zeros(5),method="SLSQP",constraints=[dict(type="ineq",fun=lambda y:1-y@y)],
                         options=dict(maxiter=150,ftol=1e-11))
            y=res.x/max(1,np.linalg.norm(res.x))
            templates.append((Q@(y*(1-1e-7)),sign))
    out=[]
    for row in data["rows"]:
        if row["passed"] or row["repeat"]:continue
        eps=F(row["epsilon"]);c=Check(eps);k=Kernel(eps)
        candidates=[]
        for z,sign in templates:
            xr=[F(f"{v:.14f}") for v in z]
            er=[s*eps/scale for s,scale in zip(sign,[1,max(F("38.1"),F.from_float(38.1)),max(F("38.1"),F.from_float(38.1))])]
            gap=gaps(z,[0,0,0])-gaps(z,list(map(float,er)))
            candidates.append((gap,xr+er))
        witnesses=[]
        for side in [0,1]:
            _,coords=max(candidates,key=lambda p:p[0][side])
            exact=[q(v) for v in coords];fixed=[pt(v) for v in coords]
            dom=c.domain(exact);fxdom=k.domain(fixed)
            assert dom<c.r and fxdom.h<k.r.l
            assert abs(coords[5])<=eps and max(map(abs,coords[6:]))*max(F("38.1"),F.from_float(38.1))<=eps
            dl=c.direct(exact)[side];fl=k.losses(fixed)[side]
            assert fraction(dl.a._mpi_[0])<=F(fl.h,S) and fraction(dl.b._mpi_[1])>=F(fl.l,S)
            if dl.a>10 and fl.l>10*S:
                witnesses.append(dict(side=["lower","upper"][side],coords=[str(x) for x in coords],
                    domain_ratio_upper=float(np.nextafter(float((dom/c.r).b),np.inf)),
                    loss_lower_N=float(np.nextafter(float(dl.a),-np.inf)),loss_upper_N=float(np.nextafter(float(dl.b),np.inf)),
                    fixed_lower_N=str(F(fl.l,S))))
        out.append(dict(epsilon=str(eps),witnesses=witnesses))
    writej(H/"counterexamples.json",out)
    for r in out:print(r["epsilon"],[(v["side"],round(v["loss_lower_N"],6)) for v in r["witnesses"]])
if __name__=="__main__":run()
