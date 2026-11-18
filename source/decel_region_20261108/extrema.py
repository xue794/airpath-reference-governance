"""Numerical proposal only. Every strict negative is checked separately."""
from model import *
from scipy.optimize import minimize
from numpy.polynomial.polynomial import polyval,polyder
def margins(t,e,D,uniform=False):
    v,n,a=center(t,D);dv,ga,dn,ng,nc=e;v+=dv;n+=dn
    u=a-R['kVp']*dv-R['kVd']*dn;ug=-R['kGammaP']*ga-R['kGammaD']*ng;uc=-R['kChiD']*nc
    s,c=np.sin(ga),np.cos(ga);ag=v*ng+g*c;ac=v*c*nc;rr=np.hypot(ag,ac)
    dag=n*ng+v*ug-g*s*ng;dac=(n*c-v*s*ng)*nc+v*c*uc
    X=m*(n+g*s);Xd=m*(u+g*c*ng);N=m*rr;Nd=m*(ag*dag+ac*dac)/rr
    def aero(a):
        q=K*v*v;Dg=q*polyval(a,cd);L=q*polyval(a,cl);Da=q*polyval(a,polyder(cd));ta=np.tan(a)
        return L+(Dg+X)*ta, q*polyval(a,polyder(cl))+Da*ta+(Dg+X)/np.cos(a)**2,2*(L+Dg*ta)/v,ta,Dg,Da
    def residual(a):return K*v*v*(polyval(a,cl)+polyval(a,cd)*np.tan(a))+X*np.tan(a)-N
    al=brentq(residual,-.05,.23);M,Ma,Mv,ta,drag,Da=aero(al)
    adot=(Nd-Mv*n-ta*Xd)/Ma;T=(drag+X)/np.cos(al)
    td=(2*drag/v*n+Da*adot+Xd)/np.cos(al)+T*ta*adot
    M,_,Mv,ta,_,_=aero(R['alphaHi']);hp=M-N-10
    return [T-3,hp,T-3+td,hp+Mv*n+ta*Xd-Nd]
def search(D,r):
    L=np.linalg.cholesky(P).T;M=np.linalg.inv(L)*np.sqrt(r)*.99999
    result=[]
    for j,start in [(0,3+D/2),(1,3+D),(2,3+D/2),(3,3+D)]:
        def fun(y):return margins(y[0],M@y[1:],D)[j]
        res=minimize(fun,np.r_[start,np.zeros(5)],method='SLSQP',bounds=[(3,35)]+[(-1,1)]*5,
            constraints={'type':'ineq','fun':lambda y:1-y[1:]@y[1:]},options={'ftol':1e-10,'maxiter':120})
        t=res.x[0];e=M@res.x[1:];result.append(dict(t=t,e=e.tolist(),VB=e@P@e,item=['hTmin','hNplus','phiTmin','phiNplus'][j],value=fun(res.x)))
    return result
if __name__=='__main__':
    for r in [.01,.003]:
        a=search(22,r);save(OUT/f'extrema_{r}.json',a);print(r,[(s['item'],s['value']) for s in a])
