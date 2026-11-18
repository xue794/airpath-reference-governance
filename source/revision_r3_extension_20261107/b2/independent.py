"""Independent 65-decimal interval replay from a heap path and rational endpoints."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sys
from boxes import H,MODEL
sys.path.insert(0,str(H.parents[1]/"revision_wp3_20261107/vendor"))
import mpmath as mp
iv=mp.iv;iv.dps=65
def q(x):
    x=F(x);return iv.mpf(x.numerator)/iv.mpf(x.denominator)
def hull(a,b):return iv.mpf([min(a.a,b.a),max(a.b,b.b)])
def poly(cs,x):
    ans=cs[-1]
    for c in cs[-2::-1]:ans=c+x*ans
    return ans
def fraction(x):
    sign,m,e,_=x
    return (-1 if sign else 1)*F(m)*F(2)**e
def outer(x):
    return [str(fraction(u)) for u in x._mpi_]

class Check:
    def __init__(self,eps,sub=False):
        b=MODEL["binary_model"]
        def number(s):return q(F.from_float(float.fromhex(s)))
        self.p=[[number(t) for t in row] for row in b["PB"]]
        self.r=number(b["rhoB"])/(4 if sub else 1)
        self.eps=F(eps)
        self.width=[F(x)/(2 if sub else 1) for x in MODEL["analytic"]["projection"]]+[self.eps,self.eps/F("38.1"),self.eps/F("38.1")]
        self.active=8 if self.eps else 5
        self.v0=number(b["Vb"]);self.m=number(b["mass"]);self.g=number(b["g"])
        self.kv,self.kg,self.kc=map(number,[b["kV"],b["kGamma"],b["kChi"]])
        self.K=number(b["density"])*number(b["S"])/2
        self.a=number(b["alphaHi"]);self.ta=iv.tan(self.a)
        cl=list(map(number,b["CL"]));cd=list(map(number,b["CD"]))
        self.bl=cl[0]; self.bu=poly(cl,self.a)+poly(cd,self.a)*self.ta
        self.caps=[q(self.eps)*q(MODEL["analytic"][k]) for k in ("Llower","Lupper")]

    def box(self,path):
        bits=bin(path)[3:]; bounds=[]
        for j,w in enumerate(self.width):
            seq=bits[j::self.active] if j<self.active else ""
            n=1<<len(seq);k=int(seq,2) if seq else 0
            bounds.append(((2*k-n)*w/n,(2*k+2-n)*w/n))
        return bounds
    def bounds(self,p):return [iv.mpf([q(l).a,q(h).b]) for l,h in self.box(p)]
    def domain(self,z):
        terms=self.p[4][4]*z[4]**2
        for i,j in [(0,2),(1,3)]:
            a,b,c=self.p[i][i],self.p[i][j],self.p[j][j]
            terms+=a*(z[i]+b/a*z[j])**2+(c-b**2/a)*z[j]**2
        return terms
    def coeff_losses(self,z):
        d,g,nv,ng,nc,x,y,c=z
        vf=self.v0+d;v=vf+x;ga=g+y
        co=iv.cos(ga);cf=iv.cos(g)
        ag0=vf*ng+self.g*cf;ac0=vf*cf*nc
        ag1=v*(ng-self.kg*y)+self.g*co;ac1=v*co*(nc-self.kc*c)
        assert ag0>0 and ag1>0
        slope=hull(ac0/ag0,ac1/ag1)
        axial=1/iv.sqrt(1+slope**2)
        f=lambda t:t/iv.sqrt(1+t**2)
        lateral=iv.mpf([f(slope.a).a,f(slope.b).b])
        arc=hull(g,ga)
        dg=-iv.sin(arc);ds=iv.cos(arc)
        N_v=self.m*(axial*ng+lateral*nc*co)
        N_g=self.m*(dg*(axial*self.g+lateral*vf*nc)-axial*v*self.kg)
        # Preserve the same distribution when comparing with stored enclosures.
        ng_direct=self.m*(axial*(self.g*dg-v*self.kg)+lateral*nc*vf*dg)
        N_g=iv.mpf([max(N_g.a,ng_direct.a),min(N_g.b,ng_direct.b)])
        N_c=-self.m*lateral*v*co*self.kc
        dq=self.K*(2*vf+x)
        lo=(dq*self.bl-N_v)*x-N_g*y-N_c*c
        hi=(N_v-dq*self.bu+self.m*self.kv*self.ta)*x+(N_g-self.m*self.g*ds*self.ta)*y+N_c*c
        return lo,hi
    def direct(self,z):
        d,g,nv,ng,nc,x,y,c=z
        def gaps(v,ga,rv,rg,rc):
            X=self.m*(rv+self.g*iv.sin(ga))
            N=self.m*iv.sqrt((v*rg+self.g*iv.cos(ga))**2+(v*iv.cos(ga)*rc)**2)
            return N-self.K*v**2*self.bl,self.K*v**2*self.bu+X*self.ta-N
        vf=self.v0+d
        a=gaps(vf,g,nv,ng,nc);b=gaps(vf+x,g+y,nv-self.kv*x,ng-self.kg*y,nc-self.kc*c)
        return [u-v for u,v in zip(a,b)]
