"""Directed fixed-point arithmetic and the regional endpoint-loss kernel."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json

H=Path(__file__).resolve().parent
B1=H.parent/"b1"/"b1_verified.json"
MODEL=json.loads(B1.read_bytes())
S=10**56
def ceil(a,b):return -((-a)//b)
class I:
    __slots__=("l","h")
    def __init__(self,l,h):self.l=l;self.h=h
    def __add__(a,b):return I(a.l+b.l,a.h+b.h)
    def __neg__(a):return I(-a.h,-a.l)
    def __sub__(a,b):return a+(-b)
    def __mul__(a,b):
        t=(a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h)
        return I(min(t)//S,ceil(max(t),S))
    def __truediv__(a,b):
        assert b.l>0
        return a*I(S*S//b.h,ceil(S*S,b.l))
    def sq(a):
        l,h=a.l,a.h
        return I((0 if l<=0<=h else min(l*l,h*h))//S,ceil(max(l*l,h*h),S))
    def sqrt(a):
        assert a.l>=0
        return I(isqrt(a.l*S),isqrt(a.h*S)+int(isqrt(a.h*S)**2!=a.h*S))
    def hull(a,b):return I(min(a.l,b.l),max(a.h,b.h))
def pt(x):
    f=F(x);return I(f.numerator*S//f.denominator,ceil(f.numerator*S,f.denominator))
def binary(s):return pt(F.from_float(float.fromhex(s)))
ONE=pt(1); ZERO=pt(0); TWO=pt(2)

@lru_cache(8192)
def trig_point(z,cosine):
    # Consecutive partial sums of the alternating Taylor series enclose the value.
    x=F(abs(z),S);assert x<F(1,4)
    n=0 if cosine else 1
    term=F(1) if cosine else x
    value=term
    while True:
        term=-term*x*x/((n+1)*(n+2)); nxt=value+term;n+=2
        if abs(term)<F(1,10**60):break
        value=nxt
    r=pt(value).hull(pt(nxt))
    return -r if z<0 and not cosine else r
def sin(x):return I(trig_point(x.l,False).l,trig_point(x.h,False).h)
def cos(x):
    low=0 if x.l<=0<=x.h else min(abs(x.l),abs(x.h))
    return I(trig_point(max(abs(x.l),abs(x.h)),True).l,trig_point(low,True).h)
def polynomial(cs,x):
    r=cs[-1]
    for c in reversed(cs[:-1]):r=c+r*x
    return r

def direction(q):
    # q/sqrt(1+q²) increases on the real line. The other component is positive.
    def lateral(z):
        a=I(z,z);return a/(ONE+a.sq()).sqrt()
    return ONE/(ONE+q.sq()).sqrt(),I(lateral(q.l).l,lateral(q.h).h)

class Kernel:
    def __init__(self,eps,sub=False):
        b=MODEL["binary_model"];a=MODEL["analytic"]
        self.eps=F(eps);r=binary(b["rhoB"])
        self.r=r*pt(F(1,4)) if sub else r
        self.P=[[binary(s) for s in row] for row in b["PB"]]
        self.half=[pt(F(s)/(2 if sub else 1)).h for s in a["projection"]]
        self.half.extend([pt(self.eps).h]+[(pt(self.eps)/pt("38.1")).h]*2)
        self.active=5 if self.eps==0 else 8
        self.v0=binary(b["Vb"]);self.m=binary(b["mass"]);self.g=binary(b["g"])
        self.kv,self.kg,self.kc=(binary(b[k]) for k in ["kV","kGamma","kChi"])
        self.K=binary(b["density"])*binary(b["S"])/TWO
        a9=binary(b["alphaHi"]);self.tan=sin(a9)/cos(a9)
        cl=[binary(s) for s in b["CL"]];cd=[binary(s) for s in b["CD"]]
        self.BL=cl[0];self.BU=polynomial(cl,a9)+polynomial(cd,a9)*self.tan
        self.block=[]
        for j,k in [(0,2),(1,3)]:
            p,q,s=self.P[j][j],self.P[j][k],self.P[k][k]
            self.block.append((j,k,p,q/p,s-q.sq()/p))
        self.caps=[(pt(self.eps)*pt(a[key])).h for key in ["Llower","Lupper"]]

    @lru_cache(8192)
    def axis(self,j,index,depth):
        n=2**depth;h=self.half[j]
        return I((2*index-n)*h//n,ceil((2*index+2-n)*h,n))

    def box(self,path):
        depth=path.bit_length()-1;bits=path-(1<<depth)
        indices=[0]*8;levels=[0]*8
        for k in range(depth):
            j=k%self.active
            indices[j]=2*indices[j]+((bits>>(depth-1-k))&1);levels[j]+=1
        return [self.axis(j,indices[j],levels[j]) for j in range(8)]

    def domain(self,z):
        result=self.P[4][4]*z[4].sq()
        for j,k,p,qp,schur in self.block:
            result=result+p*(z[j]+qp*z[k]).sq()+schur*z[k].sq()
        return result

    def losses(self,z):
        dv,gf,nv,ng,nc,ev,eg,ec=z
        vf=self.v0+dv;v=vf+ev;gamma=gf+eg
        cg=cos(gamma);cf=cos(gf)
        agf=vf*ng+self.g*cf;acf=vf*cf*nc
        ag=v*(ng-self.kg*eg)+self.g*cg;ac=v*cg*(nc-self.kc*ec)
        assert ag.l>0 and agf.l>0
        wa,wh=direction((acf/agf).hull(ac/ag))
        mid_gamma=gf.hull(gamma)
        cdiff=-sin(mid_gamma);sdiff=cos(mid_gamma)
        # Delta N is a weighted sum of component differences. For each error
        # coordinate combine the Delta N and endpoint-map coefficients first.
        dvN=self.m*(wa*ng+wh*nc*cg)
        dgN=self.m*(wa*(self.g*cdiff-v*self.kg)+wh*nc*vf*cdiff)
        dcN=-(self.m*wh*v*cg*self.kc)
        dq=self.K*(TWO*vf+ev)
        lower=(dq*self.BL-dvN)*ev-dgN*eg-dcN*ec
        upper=(dvN-dq*self.BU+self.m*self.kv*self.tan)*ev+(dgN-self.m*self.g*sdiff*self.tan)*eg+dcN*ec
        return lower,upper

    def node(self,path):
        z=self.box(path);dom=self.domain(z)
        if dom.l>self.r.h:return "P",0,0,dom.l
        losses=self.losses(z)
        u=[min(cap,max(0,v.h)) for cap,v in zip(self.caps,losses)]
        return ("A" if max(u)<10*S else "U"),*u,dom.l
