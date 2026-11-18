import numpy as np
from math import sin, cos, atan2, asin, radians, degrees
from scipy.integrate import solve_ivp
from scipy.optimize import root

ft=0.3048; ft2=ft*ft; lbf=4.4482216152605; slug=14.59390293720636; slugft2=slug*ft2
g=9.80665; W=49.6*lbf; m=W/g; S=5.902*ft2; cbar=0.915*ft; b=6.849*ft
Ixx=1.327*slugft2; Iyy=4.254*slugft2; Izz=5.454*slugft2; Ixz=0.120*slugft2
I=np.array([[Ixx,0,-Ixz],[0,Iyy,0],[-Ixz,0,Izz]],float)
th=np.array([0.019,-0.078,-27.42,0.293,3.420,288.2,-0.040,1.819,-355.3,-6.563,-1.003,0.033,0.952,-0.009,0.253,0.016,5.343,30.78,0.396,12.03,0.506,-36.30,46.13,-0.109,-0.366,0.061,-0.079,0.021,0.182,-1.782,-44.34,-1.785,374.0,-1748.0,2.439,1949.0,-0.038,0.803,0.183,-0.022,-0.405,-0.009,-0.129,0.184,-0.377])

def rho_isa(h):
    h=max(0.0,h); T0=288.15; p0=101325.; L=.0065; R=287.05287
    T=T0-L*h; p=p0*(T/T0)**(g/(R*L)); return p/(R*T)

def coef(alpha,beta,om,d,V):
    p,q,r=om; da,de,dr=d; pt=b*p/(2*V); qt=cbar*q/(2*V); rt=b*r/(2*V); t=th
    CD=t[0]+t[1]*alpha+t[2]*alpha*qt+t[3]*alpha*de+t[4]*alpha**2+t[5]*alpha**2*qt+t[6]*alpha**2*de+t[7]*alpha**3+t[8]*alpha**3*qt+t[9]*alpha**4
    CY=t[10]*beta+t[11]*pt+t[12]*rt+t[13]*da+t[14]*dr
    CL=t[15]+t[16]*alpha+t[17]*qt+t[18]*de+t[19]*alpha*qt+t[20]*alpha**2+t[21]*alpha**3+t[22]*alpha**4
    Cl=t[23]*beta+t[24]*pt+t[25]*rt+t[26]*da+t[27]*dr
    Cm=t[28]+t[29]*alpha+t[30]*qt+t[31]*de+t[32]*alpha*qt+t[33]*alpha**2*qt+t[34]*alpha**2*de+t[35]*alpha**3*qt+t[36]*alpha**3*de+t[37]*alpha**4
    Cn=t[38]*beta+t[39]*pt+t[40]*rt+t[41]*da+t[42]*dr+t[43]*beta**2+t[44]*beta**3
    return np.array([CD,CY,CL,Cl,Cm,Cn])

def euler_dcm(phi,theta,psi):
    cp,sp=cos(phi),sin(phi); ct,st=cos(theta),sin(theta); cy,sy=cos(psi),sin(psi)
    return np.array([[cy*ct,cy*st*sp-sy*cp,cy*st*cp+sy*sp],[sy*ct,sy*st*sp+cy*cp,sy*st*cp-cy*sp],[-st,ct*sp,ct*cp]])

def euler_quat(phi,theta,psi):
    hp=.5*phi; ht=.5*theta; hs=.5*psi
    cp,sp=cos(hp),sin(hp); ct,st=cos(ht),sin(ht); cs,ss=cos(hs),sin(hs)
    q=np.array([cs*ct*cp+ss*st*sp, cs*ct*sp-ss*st*cp, cs*st*cp+ss*ct*sp, ss*ct*cp-cs*st*sp])
    return q/np.linalg.norm(q)

def quat_dcm(q):
    q=q/np.linalg.norm(q); w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)], [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)], [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])

def quat_dot(q,om):
    w=q[0]; v=q[1:]; return .5*np.r_[-v@om,w*om+np.cross(v,om)]

def aero(vg,om,d,C,wI,h):
    va=vg-C.T@wI; V=np.linalg.norm(va); a=atan2(va[2],va[0]); be=asin(np.clip(va[1]/V,-1,1)); cc=coef(a,be,om,d,V); qb=.5*rho_isa(h)*V*V
    CD,CY,CL,Cl,Cm,Cn=cc
    f=qb*S*np.array([-CD*cos(a)+CL*sin(a),CY,-CD*sin(a)-CL*cos(a)])
    tau=qb*S*np.array([b*Cl,cbar*Cm,b*Cn])
    return f,tau,va

h=1200*ft; V=125*ft

def trimres(y):
    a,de,T=y; vg=V*np.array([cos(a),0,sin(a)]); C=euler_dcm(0,a,0); f,tau,_=aero(vg,np.zeros(3),np.array([0,de,0]),C,np.zeros(3),h); z=C.T@np.array([0,0,1.]); acc=g*z+(f+np.array([T,0,0]))/m
    return np.array([acc[0],acc[2],tau[1]/(m*cbar)])
tr=root(trimres,[radians(5),radians(1),18.5]).x

def wind0(t): return np.zeros(3),np.zeros(3)
def windconst(t): return np.array([8.,-4.,0.]),np.zeros(3)
def windtv(t):
    w=np.array([4+2*np.sin(.4*t),-3+1.5*np.cos(.7*t),.5*np.sin(.3*t)])
    wd=np.array([.8*np.cos(.4*t),-1.05*np.sin(.7*t),.15*np.cos(.3*t)])
    return w,wd

delta=np.array([0,tr[1],0]); T=tr[2]
def dynE(t,x,wind=wind0):
    vg=x[3:6]; phi,theta,psi=x[6:9]; om=x[9:12]; C=euler_dcm(phi,theta,psi); w,_=wind(t); f,tau,_=aero(vg,om,delta,C,w,-x[2]); z=C.T@np.array([0,0,1.]); vgd=-np.cross(om,vg)+g*z+(f+np.array([T,0,0]))/m; od=np.linalg.solve(I,tau-np.cross(om,I@om)); p,q,r=om; ed=np.array([p+q*sin(phi)*np.tan(theta)+r*cos(phi)*np.tan(theta), q*cos(phi)-r*sin(phi),(q*sin(phi)+r*cos(phi))/cos(theta)]); return np.r_[C@vg,vgd,ed,od]
def dynQ(t,x,wind=wind0,delt=delta,Tin=T):
    vg=x[3:6]; q=x[6:10]; om=x[10:13]; C=quat_dcm(q); w,_=wind(t); f,tau,_=aero(vg,om,delt,C,w,-x[2]); z=C.T@np.array([0,0,1.]); vgd=-np.cross(om,vg)+g*z+(f+np.array([Tin,0,0]))/m; od=np.linalg.solve(I,tau-np.cross(om,I@om)); return np.r_[C@vg,vgd,quat_dot(q,om),od]
def dcm_dist(A,B): return np.arccos(np.clip((np.trace(A.T@B)-1)/2,-1,1))

# P4
phi0=radians(2); theta0=tr[0]+radians(1); psi0=0; vg0=np.array([V*cos(tr[0]),0,V*sin(tr[0])]); q0=euler_quat(phi0,theta0,psi0)
xE0=np.r_[0,0,-h,vg0,phi0,theta0,psi0,0,0,0]; xQ0=np.r_[0,0,-h,vg0,q0,0,0,0]
tg=np.linspace(0,40,801)
SE=solve_ivp(lambda t,x:dynE(t,x),[0,40],xE0,t_eval=tg,rtol=1e-10,atol=1e-11)
SQ=solve_ivp(lambda t,x:dynQ(t,x),[0,40],xQ0,t_eval=tg,rtol=1e-10,atol=1e-11)
errs=[0,0,0,0,0]
for k,t in enumerate(tg):
    xe=SE.y[:,k]; xq=SQ.y[:,k]; CE=euler_dcm(*xe[6:9]); CQ=quat_dcm(xq[6:10]);
    errs[0]=max(errs[0],np.linalg.norm(xe[:3]-xq[:3])); errs[1]=max(errs[1],np.linalg.norm(xe[3:6]-xq[3:6])); errs[2]=max(errs[2],np.linalg.norm(xe[9:12]-xq[10:13])); errs[3]=max(errs[3],dcm_dist(CE,CQ)); errs[4]=max(errs[4],abs(np.linalg.norm(xq[6:10])-1))
print('P4 pos %.3e vel %.3e rate %.3e att %.3e qnorm %.3e'%tuple(errs))

# P5
C0=quat_dcm(q0); w0=np.array([8.,-4.,0.]); xQa=xQ0.copy(); xQw=xQ0.copy(); xQw[3:6]=vg0+C0.T@w0; tg2=np.linspace(0,25,501)
S0=solve_ivp(lambda t,x:dynQ(t,x,wind0),[0,25],xQa,t_eval=tg2,rtol=1e-10,atol=1e-11)
SW=solve_ivp(lambda t,x:dynQ(t,x,windconst),[0,25],xQw,t_eval=tg2,rtol=1e-10,atol=1e-11)
e5=[0,0,0,0]
for k,t in enumerate(tg2):
    a=S0.y[:,k]; z=SW.y[:,k]; C_a=quat_dcm(a[6:10]); C_z=quat_dcm(z[6:10]); va_a=a[3:6]; va_z=z[3:6]-C_z.T@w0
    e5[0]=max(e5[0],np.linalg.norm(va_a-va_z)); e5[1]=max(e5[1],dcm_dist(C_a,C_z)); e5[2]=max(e5[2],np.linalg.norm(a[10:13]-z[10:13])); e5[3]=max(e5[3],np.linalg.norm(z[:3]-a[:3]-w0*t))
print('P5 air %.3e att %.3e rate %.3e gal %.3e'%tuple(e5))

# P6
tstar=3.7; qstar=euler_quat(radians(8),radians(4),radians(20)); Cstar=quat_dcm(qstar); wstar,wdstar=windtv(tstar); vastar=np.array([35.,1.2,3.1]); vgstar=vastar+Cstar.T@wstar; omstar=np.radians([6.,-3.,4.]); dstar=np.array([radians(1.5),tr[1]+radians(.4),radians(-2)]); Tstar=T+2.; xstar=np.r_[120,-30,-350,vgstar,qstar,omstar]; xd=dynQ(tstar,xstar,windtv,dstar,Tstar); f,_,va=aero(vgstar,omstar,dstar,Cstar,wstar,350); zeta=Cstar.T@np.array([0,0,1.]); v1=xd[3:6]+np.cross(omstar,Cstar.T@wstar)-Cstar.T@wdstar; v2=-np.cross(omstar,va)+g*zeta+(f+np.array([Tstar,0,0]))/m-Cstar.T@wdstar; e6=np.linalg.norm(v1-v2)
print('P6 identity %.3e'%e6)
