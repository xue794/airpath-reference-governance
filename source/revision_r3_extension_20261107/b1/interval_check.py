"""Independent direct enclosure on the projection box using 65-decimal outward intervals."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"revision_wp3_20261107/vendor"))
import mpmath as mp

def check(C,first,old=None):
    I=mp.iv
    I.dps=65
    def point(v):
        n,d=float(v).as_integer_ratio()
        return I.mpf(n)/I.mpf(d)
    def sym(v):
        return I.mpf([-v.b,v.b])
    def hull(a,b):
        return I.mpf([min(a.a,b.a),max(a.b,b.b)])
    def merge(v,prior):
        cur=point(v)
        if old is None:return cur
        if isinstance(prior,str):prev=I.mpf(prior)
        else:prev=I.mpf([I.mpf(prior.lo).a,I.mpf(prior.hi).b])/I.mpf(old.SCALE)
        return hull(cur,prev)
    def poly(co,x):
        result=co[-1]
        for a in reversed(co[:-1]): result=result*x+a
        return result
    def norm(a,b): return I.sqrt(a**2+b**2)
    def up(a):
        from math import nextafter,inf
        return 0.0 if a.b==0 else nextafter(float(a.b),inf)
    def down(a):
        from math import nextafter,inf
        return nextafter(float(a.a),-inf)
    def require(b):
        if not b: raise ArithmeticError("B1 interval condition failed")
    R=C["R"]; P=R["P"]; A=C["A"]; M=R["PB"]
    for i in range(5):
        for j in range(5):
            if i!=j and {i,j} not in ({0,2},{1,3}): require(M[i,j]==0)
    r=merge(R["rhoB"],".10")
    expected=[None]*5
    for i,j,prev in [(0,2,old.PV if old else (None,)*3),(1,3,old.PG if old else (None,)*3)]:
        a,b,c=[merge(v,legacy) for v,legacy in zip([M[i,i],M[i,j],M[j,j]],prev)]
        require(M[i,j]==M[j,i] and a>0 and a*c-b**2>0)
        expected[i]=I.sqrt(r*c/(a*c-b**2)); expected[j]=I.sqrt(r*a/(a*c-b**2))
    require(M[4,4]>0)
    expected[4]=I.sqrt(r/merge(M[4,4],old.PC if old else None))
    if old is None:
        for w,upper in zip(expected,first["projection"]):require(w<I.mpf(upper))
    dv,gf,nv,ng,nc=[sym(v) for v in expected]
    vf=merge(R["Vb"],"38.1")+dv
    cl=[merge(P["theta"][j],s) for j,s in zip([15,16,20,21,22],[".016","5.343",".506","-36.30","46.13"])]
    cd=[merge(P["theta"][j],s) for j,s in zip([0,1,4,7,9],[".019","-.078","3.420","1.819","-6.563"])]
    m=merge(P["m"],old.m if old else None);g=merge(R["g"],"9.80665")
    K=point(R["rho"])*point(P["S"])/2
    if old:K=hull(K,I.mpf([old.C.lo,old.C.hi])/I.mpf(old.SCALE))
    kv,kg,kc=[merge(A[n],s) for n,s in zip(["kV","kGamma","kChi"],[".8","1.2","1.0"])]
    lo,hi=point(R["alphaLo"]),merge(R["alphaHi"],old.A9 if old else None); angles=I.mpf([lo.a,hi.b])
    Xf=m*(nv+g*I.sin(gf));agf=vf*ng+g*I.cos(gf);acf=vf*I.cos(gf)*nc
    nf=norm(agf,acf)
    def B(a): return poly(cl,a)+poly(cd,a)*I.tan(a)
    def fmap(a): return K*vf**2*B(a)+Xf*I.tan(a)
    gfl=m*nf-fmap(lo);gfu=fmap(hi)-m*nf
    require(gfl>10 and gfu>10)
    checks=[]
    for demo in first["examples"]:
        radius=I.mpf(demo["epsilon_mps"]);e=sym(radius)
        # Hull both exact decimal scale and the frozen binary64 scale.
        ang_scale=I.mpf([I.mpf("38.1").a,point(R["Vb"]).b])
        eg=e/ang_scale;ec=eg;v=vf+e;gam=gf+eg
        require(v>0 and abs(gam)<I.pi/2)
        d_cos=-2*I.sin(gf+eg/2)*I.sin(eg/2)
        d_sin=2*I.cos(gf+eg/2)*I.sin(eg/2)
        da=e*ng-v*kg*eg+g*d_cos
        db=(e*I.cos(gam)+vf*d_cos)*nc-v*I.cos(gam)*kc*ec
        ag=v*(ng-kg*eg)+g*I.cos(gam)
        ac=v*I.cos(gam)*(nc-kc*ec)
        den=nf+norm(ag,ac);require(den>0)
        dn_ratio=m*(da*(ag+agf)+db*(ac+acf))/den
        # Both the rationalized identity and triangle inequality bound the same quantity.
        dn_norm=sym(m*norm(da,db))
        dn=I.mpf([max(dn_norm.a,dn_ratio.a),min(dn_norm.b,dn_ratio.b)])
        dX=m*(-kv*e+g*d_sin)
        def dm(a):return K*e*(2*vf+e)*B(a)+dX*I.tan(a)
        losses=[dm(lo)-dn,dn-dm(hi)]
        X=m*(nv-kv*e+g*I.sin(gam));q=K*v**2
        Mprime=q*(poly([i*cl[i] for i in range(1,5)],angles)+
                  poly([i*cd[i] for i in range(1,5)],angles)*I.tan(angles))+(q*poly(cd,angles)+X)/I.cos(angles)**2
        require(Mprime>0)
        require(all(loss.b<=I.mpf(demo[key]).b for loss,key in zip(losses,["U_lower_N","U_upper_N"])))
        if float(radius.b)<=float(first["radius"]):
            require(all(v.b<10 for v in losses))
        if demo["epsilon_mps"]=="0.000000":require(all(v.a==0 and v.b==0 for v in losses))
        # Conservative decimal limits in the primary proof are checked independently too.
        require(gfl>I.mpf(first["raw_ref_lower"]) and gfu>I.mpf(first["raw_ref_upper"]))
        checks.append(dict(epsilon=demo["epsilon_mps"],loss_lower_upper_N=up(losses[0]),loss_upper_upper_N=up(losses[1]),
            Mprime_lower_N_per_rad=down(Mprime),reference_lower_N=down(gfl),reference_upper_N=down(gfu),
            denominator_lower_mps2=down(den),exact_zero=bool(radius==0),
            projection_halfwidths_up=[up(v) for v in expected]))
    return checks
