import numpy as np
from scipy.integrate import solve_ivp

pi=np.pi
d2r=pi/180
L=2*np.eye(3)
w0=np.array([8.,-4.,1.])
what0=np.zeros(3)
ng_bar=0.25
theta_bar=.35*d2r
nV_bar=.20
na_bar=.20*d2r
nb_bar=.15*d2r
nz_bar=np.linalg.norm([nV_bar,na_bar,nb_bar])
Vbar=45.
LT=max(1.,Vbar+nV_bar)
nw_bound_generic=ng_bar+2*Vbar*np.sin(theta_bar/2)+LT*nz_bar
nw_bound=ng_bar+2*Vbar*np.sin(theta_bar/2)+nV_bar+Vbar*(na_bar+nb_bar)

def Rz(a):
    c,s=np.cos(a),np.sin(a); return np.array([[c,-s,0],[s,c,0],[0,0,1.]])
def Ry(a):
    c,s=np.cos(a),np.sin(a); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rx(a):
    c,s=np.cos(a),np.sin(a); return np.array([[1,0,0],[0,c,-s],[0,s,c]])
def C321(phi,theta,psi): return Rz(psi)@Ry(theta)@Rx(phi)
def hat(a):
    x,y,z=a; return np.array([[0,-z,y],[z,0,-x],[-y,x,0.]])
def rotvec(rv):
    th=np.linalg.norm(rv)
    if th<1e-12:
        K=hat(rv); return np.eye(3)+K+0.5*K@K
    K=hat(rv/th); return np.eye(3)+np.sin(th)*K+(1-np.cos(th))*K@K

def airdata_to_vab(V,a,b): return V*np.array([np.cos(a)*np.cos(b),np.sin(b),np.sin(a)*np.cos(b)])
def vab_to_airdata(v):
    V=np.linalg.norm(v); return np.array([V,np.arctan2(v[2],v[0]),np.arcsin(np.clip(v[1]/V,-1,1))])
def wind(t):
    w=np.array([6+2*np.sin(.25*t), -3+1.5*np.sin(.4*t+.5), .8*np.sin(.2*t)])
    wd=np.array([.5*np.cos(.25*t), .6*np.cos(.4*t+.5), .16*np.cos(.2*t)])
    return w,wd

def noise(t):
    d=np.array([np.sin(.71*t+.2),np.cos(.43*t-.1),.6*np.sin(.31*t+.7)])
    d/=np.linalg.norm(d); ampg=ng_bar*(.65+.35*np.sin(.17*t)**2); ng=ampg*d
    a=np.array([np.cos(.29*t),np.sin(.37*t+.4),.5*np.cos(.19*t-.2)])
    a/=np.linalg.norm(a); th=theta_bar*(.55+.45*np.sin(.23*t+.1)**2); dtheta=th*a
    nz=np.array([nV_bar*np.sin(.83*t),na_bar*np.cos(.59*t+.3),nb_bar*np.sin(.47*t-.5)])
    return ng,dtheta,nz

def measurement(C,vgB,w,t,noisy):
    vgI=C@vgB; vaB=vgB-C.T@w; z=vab_to_airdata(vaB)
    if noisy:
        ng,dtheta,nz=noise(t); vgIm=vgI+ng; Cm=C@rotvec(dtheta); zm=z+nz
    else:
        ng=np.zeros(3); dtheta=np.zeros(3); nz=np.zeros(3); vgIm=vgI; Cm=C; zm=z
    yw=vgIm-Cm@airdata_to_vab(*zm)
    return yw,yw-w,ng,dtheta,nz

# P1
TRalpha=4.931803178*d2r
maxRecon=0
for k in range(101):
    t=.05*k
    C=C321(7*d2r*np.sin(.41*t),TRalpha+3*d2r*np.cos(.27*t),25*d2r*np.sin(.19*t))
    w,_=wind(t); va=np.array([34+2*np.sin(.3*t),1.2*np.sin(.5*t),3+.6*np.cos(.2*t)]); vg=va+C.T@w
    yw,*_=measurement(C,vg,w,t,False); maxRecon=max(maxRecon,np.linalg.norm(yw-w))

# P2
te=np.linspace(0,5,501)
sol=solve_ivp(lambda t,x:L@(w0-x),(0,5),what0,t_eval=te,rtol=1e-11,atol=1e-12)
what=sol.y.T; e0=what0-w0; errAnal=0; rel=[]
from scipy.linalg import expm
for t,x in zip(te,what):
    ea=expm(-L*t)@e0; en=x-w0; errAnal=max(errAnal,np.linalg.norm(en-ea)); rel.append(np.linalg.norm(en)/np.linalg.norm(e0))
idx=np.where(np.array(rel)<=.05)[0][0]; t95=te[idx]; t95th=-np.log(.05)/2

# P3
T=np.linspace(0,120,6001)
sol=solve_ivp(lambda t,x:L@(wind(t)[0]-x),(0,120),what0,t_eval=T,rtol=1e-11,atol=1e-12)
WH=sol.y.T; W=np.array([wind(t)[0] for t in T]); E=W-WH
idx=(T>=60)&(T<=120)
Om=np.array([.25,.4,.2]); Atr=np.array([2.,1.5,.8]); Afit=[]; Ath=[]
for i in range(3):
    X=np.column_stack([np.ones(idx.sum()),np.sin(Om[i]*T[idx]),np.cos(Om[i]*T[idx])])
    c=np.linalg.lstsq(X,E[idx,i],rcond=None)[0]; Afit.append(np.hypot(c[1],c[2])); Ath.append(Atr[i]*Om[i]/np.sqrt(4+Om[i]**2))
Afit=np.array(Afit); Ath=np.array(Ath); ampRel=np.abs(Afit-Ath)/Ath
rmse=np.sqrt(np.mean(E[idx]**2,axis=0)); linf=np.max(np.abs(E[idx]),axis=0)

# P4
maxNw=maxRatio=maxNg=maxTheta=maxNz=0
for k in range(1001):
    t=.03*k; C=C321(8*d2r*np.sin(.17*t),TRalpha+4*d2r*np.cos(.13*t),30*d2r*np.sin(.11*t)); w,_=wind(t)
    va=np.array([38+3*np.sin(.09*t),1*np.sin(.21*t),3.2+.8*np.cos(.15*t)]); vg=va+C.T@w
    yw,nw,ng,dt,nz=measurement(C,vg,w,t,True)
    maxNw=max(maxNw,np.linalg.norm(nw)); maxRatio=max(maxRatio,np.linalg.norm(nw)/nw_bound); maxNg=max(maxNg,np.linalg.norm(ng)); maxTheta=max(maxTheta,np.linalg.norm(dt)); maxNz=max(maxNz,np.linalg.norm(nz))

# P5
maxId=0
for k in range(501):
    t=.04*k; C=C321(6*d2r*np.sin(.23*t),TRalpha+2*d2r*np.cos(.29*t),20*d2r*np.sin(.17*t)); w,wd=wind(t)
    va=np.array([36+2*np.sin(.17*t),.8*np.sin(.31*t),3+.5*np.cos(.22*t)]); vg=va+C.T@w
    yw,nw,*_=measurement(C,vg,w,t,True)
    what=np.array([1.5*np.sin(.13*t),-.7*np.cos(.16*t),.4*np.sin(.19*t)])
    whatdot=L@(yw-what); tilde=w-what; direct=wd-whatdot; theory=-L@tilde+wd-L@nw
    maxId=max(maxId,np.linalg.norm(direct-theory))

out=f'''Section IV Wind Estimator Rev.0 independent Python reference validation\n==================================================================\nWE-P1 max ||y_w-w|| = {maxRecon:.3e} m/s\nWE-P2 max analytic mismatch = {errAnal:.3e} m/s\nWE-P2 t95 numerical = {t95:.3f} s, theory = {t95th:.3f} s\nWE-P3 fitted error amplitudes = {Afit[0]:.8f}, {Afit[1]:.8f}, {Afit[2]:.8f} m/s\nWE-P3 theory error amplitudes = {Ath[0]:.8f}, {Ath[1]:.8f}, {Ath[2]:.8f} m/s\nWE-P3 relative amplitude errors = {100*ampRel[0]:.4f}%, {100*ampRel[1]:.4f}%, {100*ampRel[2]:.4f}%\nWE-P3 steady-window RMSE = {rmse[0]:.6f}, {rmse[1]:.6f}, {rmse[2]:.6f} m/s\nWE-P3 steady-window Linf = {linf[0]:.6f}, {linf[1]:.6f}, {linf[2]:.6f} m/s\nWE-P4 max ||n_w|| = {maxNw:.6f} m/s\nWE-P4 componentwise certified bound = {nw_bound:.6f} m/s\nWE-P4 generic L_T||n_z|| bound = {nw_bound_generic:.6f} m/s\nWE-P4 max ratio = {maxRatio:.6f}\nWE-P4 realized max ||n_g|| = {maxNg:.6f} m/s\nWE-P4 realized max attitude error = {maxTheta/d2r:.6f} deg\nWE-P4 realized max ||n_z|| = {maxNz:.6f}\nWE-P5 max error-dynamics identity residual = {maxId:.3e} m/s^2\n'''
print(out)
open('/mnt/data/T2_GTM_GNA_WindEstimator_Rev0/validation_results_wind_estimator_rev0.txt','w').write(out)
