"""Mean-value enclosures over a box intersected with the exact ellipsoid."""
import numpy as np
from intervals import Box,box,trig,poly,down,up
from kernel import guide,aero,domain,K,tm,a9,nominal
from model import R,P,m,g,cl,cd
class Dual:
    __array_priority__=1100
    def __init__(self,v,d=0):self.v=box(v);self.d=box(d)
    def __add__(self,b):
        b=dual(b);return Dual(self.v+b.v,self.d+b.d)
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,-self.d)
    def __sub__(self,b):return self+-dual(b)
    def __rsub__(self,b):return dual(b)+-self
    def __mul__(self,b):
        b=dual(b);return Dual(self.v*b.v,self.d*b.v+b.d*self.v)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=dual(b);v=self.v/b.v;return Dual(v,(self.d-v*b.d)/b.v)
    def sq(self):return Dual(self.v.sq(),2*self.v*self.d)
    def sqrt(self):
        y=self.v.sqrt();return Dual(y,self.d/(2*y))
def dual(x):return x if isinstance(x,Dual) else Dual(x)
def sc(x):
    s,c=trig(x.v);return Dual(s,c*x.d),Dual(c,-s*x.d)
def pol(c,x):
    y=dual(c[-1])
    for a in c[-2::-1]:y=y*x+a
    return y
def grad(t,ex,D,phase):
    es=[Dual(x,Box(np.eye(5)[i,:,None])) for i,x in enumerate(ex)]
    v0,n0,a0=guide(t,D,phase);dv,ga,dn,ng,nc=es
    v=dv+v0;n=dn+n0;u=-dv*R['kVp']-dn*R['kVd']+a0
    ug=-ga*R['kGammaP']-ng*R['kGammaD'];uc=-nc*R['kChiD'];s,c=sc(ga)
    ag=v*ng+c*g;ac=v*c*nc;rr=(ag.sq()+ac.sq()).sqrt()
    ad=n*ng+v*ug-s*g*ng;bd=(n*c-v*s*ng)*nc+v*c*uc
    X=(n+s*g)*m;Xd=(u+c*g*ng)*m;N=rr*m;Nd=(ag*ad+ac*bd)*m/rr
    q=v.sq()*K
    def af(alpha):
        sa,ca=sc(alpha);ta=sa/ca;drag=q*pol(cd,alpha);lift=q*pol(cl,alpha)
        da=q*pol([box(i)*v for i,v in enumerate(cd) if i],alpha)
        la=q*pol([box(i)*v for i,v in enumerate(cl) if i],alpha)
        return lift+(drag+X)*ta,la+da*ta+(drag+X)/ca.sq(),(lift+drag*ta)*2/v,ta,ca,drag,da
    res={}
    for key,a,sign in [('Nminus',0,-1),('Nplus',a9,1)]:
        M,_,Mv,Mx,*_=af(dual(a))
        res['h'+key]=(M-N)*sign-10
        res['phi'+key]=(Mv*n+Mx*Xd-Nd)*sign+res['h'+key]
    aa=Box(0,R['alphaHi'])
    for _ in range(5):
        mid=box((aa.l+aa.h)*.5)
        aa=aa.cap(mid-(aero(mid,v.v,X.v)[0]-N.v)/aero(aa,v.v,X.v)[1])
    av=dual(aa);M,ma,*_=af(av)
    av.d=(N.d-M.d)/ma.v
    _,ma,mv,mx,ca,drag,da=af(av)
    adot=(Nd-mv*n-mx*Xd)/ma;T=(drag+X)/ca;Td=(drag*2/v*n+da*adot+Xd)/ca+T*mx*adot
    res.update(hTmin=T-3,hTmax=57-T,phiTmin=T-3+Td,phiTmax=57-T-Td)
    return res
def support(a,r):
    out=a[4].sq()/P[4,4]
    for i,j in [(0,2),(1,3)]:
        p,b,c=map(box,[P[i,i],P[i,j],P[j,j]])
        out+=a[i].sq()/p+(a[j]-b/p*a[i]).sq()/(c-b.sq()/p)
    return (out*r).sqrt()
def sharpen(t,e,D,phase,r):
    direct=nominal(t,e,D,phase)
    assert np.all(direct['hNplus'].l>-10) and np.all(direct['hNminus'].l>-10) and np.all(direct['Malpha'].l>0)
    mid=[box((z.l+z.h)/2) for z in e]
    base=nominal(t,mid,D,phase)
    diff=[z-c for z,c in zip(e,mid)]
    for key,y in grad(t,e,D,phase).items():
        gg=[Box(y.d.l[i],y.d.h[i]) for i in range(5)]
        a=[box((z.l+z.h)/2) for z in gg]
        w=support(a,r)
        el=Box(-w.h,w.h)-sum(ai*ci for ai,ci in zip(a,mid))
        lin=sum(ai*di for ai,di in zip(a,diff)).cap(el)
        remainder=sum((gi-ai)*di for gi,ai,di in zip(gg,a,diff))
        val=base[key]+lin+remainder
        direct[key]=Box(np.maximum(val.l,direct[key].l),np.maximum(val.h,direct[key].h))
    return direct
