"""First preselected A2 event, checked as one mathematical point of the nominal model."""
from fractions import Fraction as Q
import sys
from pathlib import Path
from reserve import binary, decimal_bound
from interval_check import mp

OLD=Path(__file__).resolve().parents[2]/"contribution_probe_20261106/backup/original_certificate"
sys.path.insert(0,str(OLD))
import interval_certificate as legacy

def old_parameters(C):
    R=C["R"]; P=R["P"]
    pairs=[("trim_V",R["Vb"],legacy.Vb),("mass",P["m"],legacy.m),("g",R["g"],legacy.g),
           ("density",R["rho"],legacy.rho),("alpha_hi",R["alphaHi"],legacy.A9),
           ("rho_B",R["rhoB"],legacy.I(".1"))]
    for (i,j),box in zip([(0,0),(0,2),(2,2),(1,1),(1,3),(3,3),(4,4)],(*legacy.PV,*legacy.PG,legacy.PC)):
        pairs.append((f"PB_{i+1}{j+1}",R["PB"][i,j],box))
    return {key:dict(runtime_hex=float(v).hex(),old_raw_interval=[iv.lo,iv.hi],scale=legacy.SCALE,
                      contains_runtime=Q(iv.lo,legacy.SCALE)<=binary(v)<=Q(iv.hi,legacy.SCALE))
            for key,v,iv in pairs}

def event(S,k):
    C=S["C"];R=C["R"];A=C["A"];P=R["P"]
    xf=list(map(binary,S["x"][k,19:25])); vf,gf,cf,nv,ng,nc=xf
    y=[binary(S["Y"][col][k]) for col in ["V","gamma","chi"]]
    ev,eg,ec=[yi-xi for yi,xi in zip(y,xf)]
    assert abs(ec)<1  # ordinary difference lies strictly within the principal arc here
    xi=[vf-binary(R["Vb"]),gf-binary(R["gammaB"]),nv,ng,nc]
    membership=sum(xi[i]*binary(R["PB"][i,j])*xi[j] for i in range(5) for j in range(5))/binary(R["rhoB"])
    width=max(abs(ev),max(Q("38.1"),binary(38.1))*max(abs(eg),abs(ec)))
    assert membership<1 and width<Q(".337374")
    iv=mp.iv; iv.dps=65
    def precise(q):return iv.mpf(q.numerator)/iv.mpf(q.denominator)
    def evaluate(make,sin,cos,tan,square,sqrt):
        # Evaluate both raw upper gaps separately at the recorded binary coordinates.
        a=make(binary(R["alphaHi"]));m=make(binary(P["m"]));g=make(binary(R["g"]))
        K=make(binary(R["rho"])*binary(P["S"])/2)
        def coeff(ids):
            h=make(Q(0))
            for j in reversed(ids):h=h*a+make(binary(P["theta"][j]))
            return h
        cl=coeff([15,16,20,21,22]);cd=coeff([0,1,4,7,9]);ta=tan(a)
        def upper(v,gam,vterm,gterm,cterm):
            X=m*(vterm+g*sin(gam)); ag=v*gterm+g*cos(gam); ac=v*cos(gam)*cterm
            N=m*sqrt(square(ag)+square(ac))
            return K*square(v)*cl+(K*square(v)*cd+X)*ta-N
        rates=[nv,ng,nc]
        feedback=[r-binary(A[n])*e for r,n,e in zip(rates,["kV","kGamma","kChi"],[ev,eg,ec])]
        f=upper(*map(make,[vf,gf,*rates])); c=upper(*map(make,[*y[:2],*feedback]))
        return f,c,f-c
    direct=evaluate(precise,iv.sin,iv.cos,iv.tan,lambda x:x**2,iv.sqrt)
    fixed=evaluate(legacy.I,legacy.sin,legacy.cos,legacy.tan,lambda x:x.sq(),lambda x:x.sqrt())
    results={}
    for name,a,b in zip(["ref_upper_N","feedback_upper_N","upper_loss_N"],direct,fixed):
        lower,upper=Q(b.lo,legacy.SCALE),Q(b.hi,legacy.SCALE)
        assert a.a>=precise(lower).a and a.b<=precise(upper).b
        results[name]=[decimal_bound(lower,10,False),decimal_bound(upper,10,True)]
    assert direct[2].a>10 and fixed[2].lo>10*legacy.SCALE and direct[1].a>0
    return dict(id=22,index_matlab=k+1,time_s=float(S["t"][k]),membership_upper=decimal_bound(membership,12,True),
        error_radius_upper_mps=decimal_bound(width,12,True),box_for_witness_mps="0.337374",**results)
