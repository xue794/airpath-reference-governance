import numpy as np
from math import sin, cos, atan2, asin, radians, degrees, sqrt
from scipy.optimize import root

ft=0.3048; ft2=ft*ft; lbf=4.4482216152605; slug=14.59390293720636; slugft2=slug*ft2
g=9.80665; W=49.6*lbf; m=W/g; S=5.902*ft2; cbar=0.915*ft; b=6.849*ft
Ixx=1.327*slugft2; Iyy=4.254*slugft2; Izz=5.454*slugft2; Ixz=0.120*slugft2
I=np.array([[Ixx,0,-Ixz],[0,Iyy,0],[-Ixz,0,Izz]],float)
th=np.array([0.019,-0.078,-27.42,0.293,3.420,288.2,-0.040,1.819,-355.3,-6.563,-1.003,0.033,0.952,-0.009,0.253,0.016,5.343,30.78,0.396,12.03,0.506,-36.30,46.13,-0.109,-0.366,0.061,-0.079,0.021,0.182,-1.782,-44.34,-1.785,374.0,-1748.0,2.439,1949.0,-0.038,0.803,0.183,-0.022,-0.405,-0.009,-0.129,0.184,-0.377])

def rho_isa(h):
    T0=288.15; p0=101325; L=0.0065; R=287.05287
    T=T0-L*h; p=p0*(T/T0)**(g/(R*L)); return p/(R*T)

def coef(alpha,beta,om,delta,V):
    p,q,r=om; da,de,dr=delta; pt=b*p/(2*V); qt=cbar*q/(2*V); rt=b*r/(2*V); t=th
    CD=t[0]+t[1]*alpha+t[2]*alpha*qt+t[3]*alpha*de+t[4]*alpha**2+t[5]*alpha**2*qt+t[6]*alpha**2*de+t[7]*alpha**3+t[8]*alpha**3*qt+t[9]*alpha**4
    CY=t[10]*beta+t[11]*pt+t[12]*rt+t[13]*da+t[14]*dr
    CL=t[15]+t[16]*alpha+t[17]*qt+t[18]*de+t[19]*alpha*qt+t[20]*alpha**2+t[21]*alpha**3+t[22]*alpha**4
    Cl=t[23]*beta+t[24]*pt+t[25]*rt+t[26]*da+t[27]*dr
    Cm=t[28]+t[29]*alpha+t[30]*qt+t[31]*de+t[32]*alpha*qt+t[33]*alpha**2*qt+t[34]*alpha**2*de+t[35]*alpha**3*qt+t[36]*alpha**3*de+t[37]*alpha**4
    Cn=t[38]*beta+t[39]*pt+t[40]*rt+t[41]*da+t[42]*dr+t[43]*beta**2+t[44]*beta**3
    return np.array([CD,CY,CL,Cl,Cm,Cn])

def dcm(phi,theta,psi):
    cp,sp=cos(phi),sin(phi); ct,st=cos(theta),sin(theta); cy,sy=cos(psi),sin(psi)
    return np.array([[cy*ct,cy*st*sp-sy*cp,cy*st*cp+sy*sp],[sy*ct,sy*st*sp+cy*cp,sy*st*cp-cy*sp],[-st,ct*sp,ct*cp]])

def aero(v,om,delta,C,h):
    V=np.linalg.norm(v); a=atan2(v[2],v[0]); be=asin(np.clip(v[1]/V,-1,1)); cc=coef(a,be,om,delta,V); qb=.5*rho_isa(h)*V*V
    CD,CY,CL,Cl,Cm,Cn=cc
    f=qb*S*np.array([-CD*cos(a)+CL*sin(a),CY,-CD*sin(a)-CL*cos(a)])
    tau=qb*S*np.array([b*Cl,cbar*Cm,b*Cn])
    return f,tau,cc

h=1200*ft; V=125*ft
C1=coef(radians(5),0,np.zeros(3),np.array([0,radians(1.8),0]),V)

def tres(y):
    a,de,T=y; vb=V*np.array([cos(a),0,sin(a)]); C=dcm(0,a,0); f,tau,_=aero(vb,np.zeros(3),np.array([0,de,0]),C,h); zeta=C.T@np.array([0,0,1.]); acc=g*zeta+(f+np.array([T,0,0]))/m
    return np.array([acc[0],acc[2],tau[1]/(m*cbar)])
tr=root(tres,np.array([radians(5),radians(1),18.5])).x

def dyn9(x):
    u,v,w,p,q,r,phi,theta,psi=x; vb=np.array([u,v,w]); om=np.array([p,q,r]); C=dcm(phi,theta,psi); f,tau,_=aero(vb,om,np.array([0,tr[1],0]),C,h); zeta=C.T@np.array([0,0,1.]); vd=-np.cross(om,vb)+g*zeta+(f+np.array([tr[2],0,0]))/m; od=np.linalg.solve(I,tau-np.cross(om,I@om)); phid=p+q*sin(phi)*np.tan(theta)+r*cos(phi)*np.tan(theta); thetad=q*cos(phi)-r*sin(phi); psid=(q*sin(phi)+r*cos(phi))/cos(theta); return np.r_[vd,od,phid,thetad,psid]
x0=np.array([V*cos(tr[0]),0,V*sin(tr[0]),0,0,0,0,tr[0],0])
A=np.zeros((9,9))
for i in range(9):
    hh=1e-6*max(1,abs(x0[i])); xp=x0.copy(); xm=x0.copy(); xp[i]+=hh; xm[i]-=hh; A[:,i]=(dyn9(xp)-dyn9(xm))/(2*hh)
ev=np.linalg.eigvals(A)
real=ev[(abs(ev.imag)<1e-8)&(abs(ev.real)>1e-8)]
roll=real[np.argmin(abs(abs(real.real)-5.44))]; spiral=real[np.argmin(abs(abs(real.real)-.0887))]
cmp=ev[ev.imag>1e-8]
def near(target): return cmp[np.argmin(abs(abs(cmp)-target))]
ph=near(.342); du=near(5.27); sh=near(6.49)
print('P1: CD=%.12f CL=%.12f Cm=%.12f'%(C1[0],C1[2],C1[4]))
print('P2: alpha=%.9f deg de=%.9f deg T=%.9f N residual=%.3e'%(degrees(tr[0]),degrees(tr[1]),tr[2],np.linalg.norm(tres(tr))))
for n,z,refw,refz in [('phugoid',ph,.342,.0452),('dutch',du,5.27,.179),('short',sh,6.49,.351)]:
    wn=abs(z); ze=-z.real/wn; print('P3 %-8s wn=%.9f ref=%.4f zeta=%.9f ref=%.4f'%(n,wn,refw,ze,refz))
print('P3 spiral=%.9f ref=.0887'%abs(spiral))
print('P3 roll=%.9f ref=5.44'%abs(roll))
