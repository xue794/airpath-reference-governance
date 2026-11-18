"""Regional lower bounds, with the moving-center backup input applied explicitly."""
import numpy as np
from intervals import Box,box,poly,trig
from model import R,A,cl,cd,m,g
K=box(.5)*R['rho']*R['P']['S']
a9=box(R['alphaHi']); s9,c9=trig(a9);tan9=s9/c9
sb,cb=trig(box(R['muMax']));tm=sb/cb
def deriv(c):return [box(a)*i for i,a in enumerate(c) if i]
def aero(a,V,X):
    sg,cg=trig(a);ta=sg/cg;q=K*V.sq();D=q*poly(cd,a);L=q*poly(cl,a);Da=q*poly(deriv(cd),a)
    N=L+(D+X)*ta;Ma=q*poly(deriv(cl),a)+Da*ta+(D+X)/cg.sq()
    return N,Ma,2*(L+D*ta)/V,ta,cg,D,Da
def guide(t,D,phase):
    a,b,h=(R['Vb'],31.25,D) if phase==1 else (31.25,R['Vb'],20)
    if phase in (0,2):return box(R['Vb'] if phase==0 else 31.25),box(0),box(0)
    q=((t-(3 if phase==1 else 35))/h).cap(Box(0,1))
    def S(z):return z*z*z*(10+z*(-15+6*z))
    s=Box(S(box(q.l)).l,S(box(q.h)).h)
    # q(1-q)=1/4-(q-1/2)^2
    w=box(.25)-(q-.5).sq()
    return a+(box(b)-a)*s,(box(b)-a)*30*w.sq()/h,(box(b)-a)*60*w*(1-2*q)/(h*h)
def domain(e):
    p=R['PB'];v=box(p[4,4])*e[4].sq()
    for i,j in [(0,2),(1,3)]:
        a,b,c=map(box,[p[i,i],p[i,j],p[j,j]])
        v+=a*(e[i]+b/a*e[j]).sq()+(c-b.sq()/a)*e[j].sq()
    return v
def nominal(t,e,D,phase):
    v0,n0,u0=guide(t,D,phase);dv,ga,dn,ng,nc=e
    v=v0+dv;nv=n0+dn;u=u0-R['kVp']*dv-R['kVd']*dn
    ug=-R['kGammaP']*ga-R['kGammaD']*ng;uc=-R['kChiD']*nc
    sg,cg=trig(ga);ag=v*ng+g*cg;ac=v*cg*nc
    agd=nv*ng+v*ug-g*sg*ng;acd=(nv*cg-v*sg*ng)*nc+v*cg*uc
    X=m*(nv+g*sg);Xd=m*(u+g*cg*ng);rr=(ag.sq()+ac.sq()).sqrt()
    N=m*rr;Nd=m*(ag*agd+ac*acd)/rr
    h=[v-R['Vmin'],R['Vmax']-v,ga+R['gammaMax'],R['gammaMax']-ga]
    rates=[nv,-nv,ng,-ng];acc=[u,-u,ug,-ug]; k=box(R['k1']); k2=box(R['k2'])
    res={}
    for j,(key,v0_) in enumerate(zip(['Vmin','Vmax','Gmin','Gmax'],h)):
        res['h'+key]=v0_;res['psi'+key]=rates[j]+k*v0_
        res['phi'+key]=acc[j]+(k+k2)*rates[j]+k*k2*v0_
    for j,(w,wd) in enumerate([(nv,u),(ng,ug),(nc,uc)]):
        for name,s in [('min',1),('max',-1)]:
            key=['nuV','nuG','nuC'][j]+name
            res['h'+key]=R['nuMax'][j]+s*w;res['phi'+key]=s*wd+R['kRate']*res['h'+key]
    for key,val,delta,rate in [('AG',ag-R['aGammaMin'],agd,R['kAGamma']),('MuPlus',tm*ag-ac,tm*agd-acd,R['kMu']),('MuMinus',tm*ag+ac,tm*agd+acd,R['kMu'])]:
        res['h'+key]=val;res['phi'+key]=delta+rate*val
    for key,a,s in [('Nminus',box(0),-1),('Nplus',a9,1)]:
        force,_,mv,mx,*_=aero(a,v,X)
        res['h'+key]=s*(force-N)-R['DeltaN']
        res['phi'+key]=s*(mv*nv+mx*Xd-Nd)+R['kN']*res['h'+key]
    aa=Box(0,R['alphaHi']);res['Malpha']=aero(aa,v,X)[1]
    for _ in range(5):
        mid=box((aa.l+aa.h)*.5)
        f=aero(mid,v,X)[0]-N;df=aero(aa,v,X)[1]
        aa=aa.cap(mid-f/df)
    assert np.all(aa.l<=aa.h)
    _,ma,mv,mx,ca,drag,da=aero(aa,v,X)
    ad=(Nd-mv*nv-mx*Xd)/ma;T=(drag+X)/ca;Td=(2*drag/v*nv+da*ad+Xd)/ca+T*mx*ad
    res.update(hTmin=T-(R['Tmin']+R['DeltaT']),hTmax=(R['Tmax']-R['DeltaT'])-T,cosAlpha=ca,forceNorm=rr)
    res['phiTmin']=Td+R['kT']*res['hTmin'];res['phiTmax']=-Td+R['kT']*res['hTmax']
    return res

def feedback(t,e,D,phase,uniform=False):
    v0,n0,_=guide(t,D,phase);dv,ga,dn,ng,nc=e;vf=v0+dv;nv=n0+dn
    # Frozen decimal angle widths and pi are enclosed, not rounded inward.
    pi=Box(np.nextafter(np.pi,-np.inf),np.nextafter(np.pi,np.inf))
    wg=box(45)*pi/18000;wc=wg if uniform else box(35)*pi/18000
    wv=box(381)/10*wg if uniform else box(.125)
    unit=Box(-1,1);ev=unit*wv;eg=unit*wg;ec=unit*wc
    v=vf+ev;gamma=ga+eg;si,co=trig(gamma);sf,cf=trig(ga)
    kv,kg,kc=map(box,[A['kV'],A['kGamma'],A['kChi']])
    af=vf*ng+g*cf;bf=vf*cf*nc;rf=(af.sq()+bf.sq()).sqrt()
    Xf=m*(nv+g*sf);X=m*(nv-kv*ev+g*si)
    ag=v*(ng-kg*eg)+g*co;ac=v*co*(nc-kc*ec);rr=(ag.sq()+ac.sq()).sqrt()
    coeff=[box(cl[0]),poly(cl,a9)+poly(cd,a9)*tan9]
    q=K*v.sq();qf=K*vf.sq()
    fb0=m*rr-q*coeff[0];fb1=q*coeff[1]+X*tan9-m*rr
    # At the upper endpoint g+ increases with eV and egamma (checked per box).
    # Its minimum over course error is attained at one of the two faces.
    dV=2*K*v*coeff[1]-m*kv*tan9-m*(ag*(ng-kg*eg)+ac*co*(nc-kc*ec))/rr
    dg=m*g*co*tan9-m*(ag*(-v*kg-g*si)-ac*v*si*(nc-kc*ec))/rr
    def losses(ev,eg,ec):
        v=vf+ev;si,co=trig(ga+eg);sh,ch=trig(ga+Box(np.minimum(0,eg.l),np.maximum(0,eg.h)))
        ag1=v*(ng-kg*eg)+g*co;ac1=v*co*(nc-kc*ec)
        da=ev*ng-v*kg*eg-g*sh*eg
        db=(ev*co-vf*sh*eg)*nc-v*co*kc*ec
        # Direction along the line segment lies between the two component ratios.
        qa=bf/af;qb=ac1/ag1
        ratios=Box(np.minimum(qa.l,qb.l),np.maximum(qa.h,qb.h))
        wa=box(1)/(1+ratios.sq()).sqrt()
        def lateral(z):return z/(1+z.sq()).sqrt()
        wh=Box(lateral(box(ratios.l)).l,lateral(box(ratios.h)).h)
        dn=m*(wa*da+wh*db);dq=K*(2*vf+ev)*ev
        return dq*coeff[0]-dn,dn-dq*coeff[1]-m*(-kv*ev+g*ch*eg)*tan9
    lower,upper=losses(ev,eg,ec)
    left=losses(-wv,-wg,-wc)[1];right=losses(-wv,-wg,wc)[1]
    bound=np.where((dV.l>0)&(dg.l>0),np.minimum(upper.h,np.maximum(left.h,right.h)),upper.h)
    upper=Box(np.minimum(upper.l,bound),bound)
    res={'payLower':10-lower,'payUpper':10-upper,'fbLo':fb0,'fbHi':fb1,
         'fbMa':aero(Box(0,R['alphaHi']),v,X)[1],'normal':ag-R['aGammaMin'],
         'bankPlus':tm*ag-ac,'bankMinus':tm*ag+ac,'fbV':v}
    aa=Box(0,R['alphaHi'])
    for _ in range(4):
        mid=box((aa.l+aa.h)*.5);aa=aa.cap(mid-(aero(mid,v,X)[0]-m*rr)/aero(aa,v,X)[1])
    assert np.all(aa.l<=aa.h)
    _,_,_,_,_,drag,_=aero(aa,v,X)
    kap=Box(trig(aa+box(7)*pi/1800)[1].l,1.)
    T=(drag+X)/kap
    res.update(kappa=kap,rawTmin=T,rawTmax=60-T)
    return res
