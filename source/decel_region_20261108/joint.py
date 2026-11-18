"""Enclose feedback axial force over the ellipsoid and the three error axes."""
import numpy as np
from centered import Dual,pol,sc,support
from kernel import guide,K,aero
from intervals import box,Box,down,up,trig
from model import A,R,cl,cd,m,g
def numerator(t,coords,D,phase):
    vc,nv,_=guide(t,D,phase)
    s=[Dual(a,Box(np.eye(8)[j,:,None])) for j,a in enumerate(coords)]
    dv,ga,dn,ng,nc,ev,eg,ec=s
    v=dv+vc+ev;sg,cg=sc(ga+eg)
    X=(dn+nv-ev*A['kV']+sg*g)*m
    ag=v*(ng-eg*A['kGamma'])+cg*g;ac=v*cg*(nc-ec*A['kChi'])
    N=(ag.sq()+ac.sq()).sqrt()*m;q=v.sq()*K
    aa=Box(0,R['alphaHi'])
    assert np.all((N.v-aero(box(0),v.v,X.v)[0]).l>0) and np.all((aero(box(R['alphaHi']),v.v,X.v)[0]-N.v).l>0)
    assert np.all(aero(aa,v.v,X.v)[1].l>0)
    for _ in range(5):
        mid=box((aa.l+aa.h)*.5)
        aa=aa.cap(mid-(aero(mid,v.v,X.v)[0]-N.v)/aero(aa,v.v,X.v)[1])
    assert np.all(aa.l<=aa.h)
    al=Dual(aa);si,co=sc(al)
    residual=q*pol(cl,al)+(q*pol(cd,al)+X)*si/co-N
    al.d=-residual.d/aero(aa,v.v,X.v)[1]
    return q*pol(cd,al)+X
def joint(t,e,D,phase,r):
    pi=Box(down(np.pi),up(np.pi));unit=Box(-1,1)
    coords=[*e,unit*.125,unit*45*pi/18000,unit*35*pi/18000]
    mid=[box((z.l+z.h)*.5) for z in coords]
    val=numerator(t,coords,D,phase);base=numerator(t,mid,D,phase).v
    gg=[Box(val.d.l[i],val.d.h[i]) for i in range(8)]
    a=[box((v.l+v.h)*.5) for v in gg]
    dx=[z-c for z,c in zip(coords,mid)];w=support(a[:5],r)
    lin=sum(ai*di for ai,di in zip(a[:5],dx[:5])).cap(Box(-w.h,w.h)-sum(ai*b for ai,b in zip(a[:5],mid[:5])))
    val=val.v.cap(base+lin+sum(ai*di for ai,di in zip(a[5:],dx[5:]))+sum((gi-ai)*di for gi,ai,di in zip(gg,a,dx)))
    # For the conditional attitude sum <=0.70 degrees, kappa>=cos(9deg+0.70deg).
    kap=Box(trig(box(R['alphaHi'])+7*pi/1800)[1].l,1.)
    T=val/kap
    return dict(rawTmin=T,rawTmax=60-T)
