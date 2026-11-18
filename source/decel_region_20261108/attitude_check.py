"""Optional raw-thrust test with zero attitude discrepancy and a legal box error."""
from model import *
sys.path.insert(0,str(OLD))
from screen_feedback import ForceCheck,guide as exact_guide,pt,I,lo,hi,wp,BOXES
def witness(D,r,typ):
    eng=ForceCheck(data['C']);mm=eng.model;w=BOXES[typ]
    # A positive speed error and negative climb error reduce axial force.
    # Choose an interior error to avoid binary/real boundary ambiguity.
    eair=np.array([float(z.mid) for z in w])*np.array([.999,-.999,.999])
    M=np.linalg.inv(np.linalg.cholesky(P).T)*np.sqrt(r)*.999
    def evaluate(y):
        t=y[0];v,n,_=center(t,D);dv,ga,dn,ng,nc=M@y[1:]
        v+=dv+eair[0];ga+=eair[1];X=m*(n+dn-A['kV']*eair[0]+g*np.sin(ga))
        N=m*np.hypot(v*(ng-A['kGamma']*eair[1])+g*np.cos(ga),v*np.cos(ga)*(nc-A['kChi']*eair[2]))
        drag=lambda a:K*v*v*np.polynomial.polynomial.polyval(a,cd)
        al=brentq(lambda a:K*v*v*np.polynomial.polynomial.polyval(a,cl)+(drag(a)+X)*np.tan(a)-N,-.06,.24)
        return (drag(al)+X)/np.cos(al)
    res=minimize(evaluate,np.r_[14.,np.zeros(5)],method='SLSQP',bounds=[(14,20)]+[(-1,1)]*5,
                 constraints={'type':'ineq','fun':lambda y:1-y[1:]@y[1:]},options={'ftol':1e-11})
    t=float(res.x[0]);off=M@res.x[1:];v0=exact_guide(t,D)
    x=[u+pt(o) for u,o in zip(v0,[*off[:2],0,*off[2:]])]
    vb=sum(pt(off[i])*pt(P[i,j])*pt(off[j]) for i in range(5) for j in range(5))
    assert hi(vb)<r and all(pt(abs(e))<z for e,z in zip(eair,w)) and 14<=t<=20
    err=list(map(pt,eair));gaps,q,X=eng.evaluate(x,err);assert min(map(lo,gaps))>0 and eng.slopes(q,X)>0
    v=x[0]+err[0];ga=x[1]+err[1];N=mm.m*I.sqrt((v*(x[4]-mm.kg*err[1])+mm.g*I.cos(ga))**2+(v*I.cos(ga)*(x[5]-mm.kc*err[2]))**2)
    a,b=mm.alo,mm.ahi
    for _ in range(49):
        mid=pt((a+b)*.5);f=q*wp.poly(mm.cl,mid)+(q*wp.poly(mm.cd,mid)+X)*I.tan(mid)-N
        if lo(f)>0:b=float(mid.mid)
        elif hi(f)<0:a=float(mid.mid)
        else:break
    al=I.mpf([a,b]);T=(q*wp.poly(mm.cd,al)+X)/I.cos(al)
    return dict(D=D,r=r,box=typ,time_s=t,offset=off.tolist(),tracking_error=eair.tolist(),VB_upper=hi(vb),
                raw_thrust=[lo(T),hi(T)],zero_attitude_discrepancy=True,counterexample=bool(hi(T)<0))
if __name__=='__main__':
    r=float(sys.argv[1]);save(OUT/f'raw_thrust_{r}.json',[witness(22,r,b) for b in BOXES])
