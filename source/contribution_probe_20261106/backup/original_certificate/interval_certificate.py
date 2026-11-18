"""Exact outward-rounded fixed-point intervals for Paper 1 reference model.
This verifies only the explicitly transcribed mathematical reference model,
not MATLAB floating-point execution or full physical closed-loop safety.
"""
from fractions import Fraction
from math import isqrt,factorial
import time,json,sys
from functools import lru_cache
SCALE=10**18

def cd(a,b):return -((-a)//b)
class I:
 __slots__=('lo','hi')
 def __init__(self,a=0,b=None,raw=False):
  if raw:self.lo,self.hi=a,b;return
  if isinstance(a,I):self.lo,self.hi=a.lo,a.hi;return
  a=Fraction(str(a));b=a if b is None else Fraction(str(b))
  self.lo=a.numerator*SCALE//a.denominator;self.hi=cd(b.numerator*SCALE,b.denominator)
 @staticmethod
 def raw(lo,hi):return I(lo,hi,True)
 def __add__(self,o):
  o=I(o);return I.raw(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I.raw(-self.hi,-self.lo)
 def __sub__(self,o):return self+-I(o)
 def __rsub__(self,o):return I(o)+-self
 def __mul__(self,o):
  o=I(o);p=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];return I.raw(min(p)//SCALE,cd(max(p),SCALE))
 __rmul__=__mul__
 def inv(self):
  if self.lo<=0<=self.hi:raise ZeroDivisionError('contains zero')
  return I.raw(SCALE*SCALE//self.hi,cd(SCALE*SCALE,self.lo))
 def __truediv__(self,o):return self*I(o).inv()
 def __rtruediv__(self,o):return I(o)*self.inv()
 def sq(self):
  a,b=self.lo,self.hi;lo=0 if a<=0<=b else min(a*a,b*b)//SCALE
  return I.raw(lo,cd(max(a*a,b*b),SCALE))
 def sqrt(self):
  assert self.lo>=0
  a=isqrt(self.lo*SCALE);b=isqrt(self.hi*SCALE)
  if b*b<self.hi*SCALE:b+=1
  return I.raw(a,b)
 def intersect(self,o):
  o=I(o);a=max(self.lo,o.lo);b=min(self.hi,o.hi)
  if a>b:raise ValueError('empty interval')
  return I.raw(a,b)
 def absmax(self):
  v=max(abs(self.lo),abs(self.hi));return I.raw(v,v)
 def midpoint(self):
  v=(self.lo+self.hi)//2;return I.raw(v,v)
 def pair(self):return [self.lo/SCALE,self.hi/SCALE]
 def __repr__(self):return str(self.pair())
 def __hash__(self):return hash((self.lo,self.hi))
 def __eq__(self,o):return isinstance(o,I) and self.lo==o.lo and self.hi==o.hi
def ipow(x,n):
 r=I(1)
 for _ in range(n):r=r*x
 return r
@lru_cache(maxsize=32768)
def sinpt(x):
 r=I(0)
 for k in range(8):r+=((-1)**k)*ipow(x,2*k+1)/factorial(2*k+1)
 e=ipow(x.absmax(),16)/factorial(16);return r+I.raw(-e.hi,e.hi)
@lru_cache(maxsize=32768)
def cospt(x):
 r=I(0)
 for k in range(8):r+=((-1)**k)*ipow(x,2*k)/factorial(2*k)
 e=ipow(x.absmax(),15)/factorial(15);return r+I.raw(-e.hi,e.hi)
def sin(x):
 assert x.absmax().hi<I('1.1').lo
 return I.raw(sinpt(I.raw(x.lo,x.lo)).lo,sinpt(I.raw(x.hi,x.hi)).hi)
def cos(x):
 assert x.absmax().hi<I('1.1').lo
 lo=cospt(x.absmax()).lo
 if x.lo<=0<=x.hi:hi=SCALE
 else:
  v=min(abs(x.lo),abs(x.hi));hi=cospt(I.raw(v,v)).hi
 return I.raw(lo,hi)
def tan(x):return sin(x)/cos(x)
def atanpt(x,n=32):
 r=I(0)
 for k in range(n):r+=((-1)**k)*ipow(x,2*k+1)/(2*k+1)
 e=ipow(x.absmax(),2*n+1)/(2*n+1);return r+I.raw(-e.hi,e.hi)
def log1(x,n=12):
 z=(x-1)/(x+1);r=I(0)
 for k in range(n):r+=2*ipow(z,2*k+1)/(2*k+1)
 az=z.absmax();e=2*ipow(az,2*n+1)/((2*n+1)*(1-az.sq()));return r+I.raw(-e.hi,e.hi)
def exps(x,n=18):
 assert x.absmax().hi<I(1).lo
 r=I(0)
 for k in range(n+1):r+=ipow(x,k)/factorial(k)
 e=3*ipow(x.absmax(),n+1)/factorial(n+1);return r+I.raw(-e.hi,e.hi)
PI=16*atanpt(I(1)/5)-4*atanpt(I(1)/239);DEG=PI/180
m=I('49.6')*I('4.4482216152605')/I('9.80665');g=I('9.80665');Vb=I('38.1')
T0=I('288.15');Tm=T0-I('.0065')*I('365.76');rho=I(101325)*exps((g/(I('287.05287')*I('.0065')))*log1(Tm/T0))/(I('287.05287')*Tm)
C=I('.5')*rho*I('5.902')*I('.3048').sq();q1=I('.04');qg=1/(5*DEG).sq();qv=I(1);qc=1/(8*DEG).sq()
def pb(a,b,q1,q2):
 a=I(a);b=I(b);p12=q1/(2*a);p22=(q2+q1/a)/(2*b);p11=b*p12+a*p22
 return p11,p12,p22
PV=pb('.6','1.4',q1,qv);PG=pb('1.2','1.8',qg,qg);PC=qc/3
# Fixed constants used repeatedly.
TM=tan(35*DEG);A9=9*DEG

def aero(a,V):
 CD=I('.019')+a*(I('-.078')+a*(I('3.420')+a*(I('1.819')-I('6.563')*a)))
 CL=I('.016')+a*(I('5.343')+a*(I('.506')+a*(I('-36.30')+I('46.13')*a)))
 CDa=I('-.078')+a*(I('6.840')+a*(I('5.457')-I('26.252')*a))
 CLa=I('5.343')+a*(I('1.012')+a*(I('-108.9')+I('184.52')*a))
 q=C*V.sq();return q*CD,q*CL,q*CDa,q*CLa

def nm(a,V,X):
 D,L,Da,La=aero(a,V);ta=tan(a);ca=cos(a)
 return L+(D+X)*ta,La+Da*ta+(D+X)/ca.sq(),2*(L+D*ta)/V,ta

def ellipse(z):
 v,ga,nv,ng,nc=z;a,b,c=PV;d,e,f=PG
 return a*(v+b/a*nv).sq()+(c-b.sq()/a)*nv.sq()+d*(ga+e/d*ng).sq()+(f-e.sq()/d)*ng.sq()+PC*nc.sq()

def evaluate(z,rootsteps=5):
 dv,ga,nv,ng,nc=z;V=Vb+dv
 uV=-I('.6')*dv-I('1.4')*nv;uG=-I('1.2')*ga-I('1.8')*ng;uC=-I('1.5')*nc
 sg=sin(ga);cg=cos(ga);ag=V*ng+g*cg;ac=V*cg*nc
 agd=nv*ng+V*uG-g*sg*ng;acd=(nv*cg-V*sg*ng)*nc+V*cg*uC
 X=m*(nv+g*sg);Xd=m*(uV+g*cg*ng);rr=(ag.sq()+ac.sq()).sqrt();N=m*rr;Nd=m*(ag*agd+ac*acd)/rr
 ML=nm(I(0),V,X);MH=nm(A9,V,X);deriv=nm(I.raw(0,A9.hi),V,X)[1]
 out={'hVmin':V-30,'hVmax':46-V,'hGmin':ga+10*DEG,'hGmax':10*DEG-ga,
 'psiVmin':nv+I('1.2')*(V-30),'psiVmax':-nv+I('1.2')*(46-V),'psiGmin':ng+I('1.2')*(ga+10*DEG),'psiGmax':-ng+I('1.2')*(10*DEG-ga),
 'hnuVmin':nv+2,'hnuVmax':2-nv,'hnuGmin':ng+8*DEG,'hnuGmax':8*DEG-ng,'hnuCmin':nc+12*DEG,'hnuCmax':12*DEG-nc,
 'hAG':ag-4,'hMuPlus':TM*ag-ac,'hMuMinus':TM*ag+ac,'hNminus':N-ML[0]-10,'hNplus':MH[0]-N-10,'N_alpha':deriv,
 'phiVmin':I('.84')*dv+nv+I('1.44')*(Vb-30),'phiVmax':-I('.84')*dv-nv+I('1.44')*(46-Vb),
 'phiGmin':I('.24')*ga+I('.6')*ng+I('1.44')*10*DEG,'phiGmax':-I('.24')*ga-I('.6')*ng+I('1.44')*10*DEG,
 'phinuVmin':-I('.6')*dv+I('.1')*nv+3,'phinuVmax':I('.6')*dv-I('.1')*nv+3,
 'phinuGmin':-I('1.2')*ga-I('.3')*ng+I('1.5')*8*DEG,'phinuGmax':I('1.2')*ga+I('.3')*ng+I('1.5')*8*DEG,
 'phinuCmin':I('1.5')*12*DEG,'phinuCmax':I('1.5')*12*DEG,'phiAG':agd+ag-4,
 'phiMuPlus':TM*(agd+ag)-(acd+ac),'phiMuMinus':TM*(agd+ag)+(acd+ac),
 'phiNminus':Nd-ML[2]*nv-ML[3]*Xd+N-ML[0]-10,'phiNplus':MH[2]*nv+MH[3]*Xd-Nd+MH[0]-N-10}
 if min(out['hNminus'].lo,out['hNplus'].lo,deriv.lo)<=0:return out,False
 aa=I.raw(0,A9.hi)
 for k in range(rootsteps):
  md=aa.midpoint();fv=nm(md,V,X)[0]-N;df=nm(aa,V,X)[1]
  if df.lo<=0:df=deriv
  new=aa.intersect(md-fv/df)
  if new.hi-new.lo >= aa.hi-aa.lo:break
  aa=new
 Ma=nm(aa,V,X);D,L,Da,La=aero(aa,V);ad=(Nd-Ma[2]*nv-Ma[3]*Xd)/Ma[1]
 T=(D+X)/cos(aa);Td=(2*D/V*nv+Da*ad+Xd)/cos(aa)+T*tan(aa)*ad
 out.update({'hTmin':T-3,'hTmax':57-T,'phiTmin':Td+T-3,'phiTmax':-Td+57-T,'cosAlpha':cos(aa)})
 return out,True

def projected(rB):
 a,b,c=PV;d,e,f=PG;det=a*c-b.sq();detg=d*f-e.sq()
 return [(rB*k).sqrt().hi for k in [c/det,f/detg,a/det,d/detg,1/PC]]

def verify(rb='0.1',maxnodes=250000):
 rB=I(rb);half=projected(rB);stack=[([I.raw(-h,h) for h in half],0)];num=0;pruned=0;accepted=0;lower={};depthmax=0;start=time.time()
 while stack:
  z,depth=stack.pop();num+=1;depthmax=max(depthmax,depth)
  if ellipse(z).lo>rB.hi:pruned+=1;continue
  try:out,valid=evaluate(z)
  except (ZeroDivisionError,ValueError):valid=False;out={}
  if valid and all(v.lo>0 for v in out.values()):
   accepted+=1
   for k,v in out.items():lower[k]=min(lower.get(k,v.lo),v.lo)
  else:
   if depth>=45 or num>=maxnodes:return {'certified':False,'rho':rb,'nodes':num,'unresolved':len(stack)+1,'depth':depth,'failed':{k:v.pair() for k,v in out.items() if v.lo<=0},'seconds':time.time()-start}
   j=max(range(5),key=lambda j:(z[j].hi-z[j].lo)/half[j]);mid=(z[j].lo+z[j].hi)//2
   z1=z.copy();z2=z.copy();z1[j]=I.raw(z[j].lo,mid);z2[j]=I.raw(mid,z[j].hi);stack.extend([(z2,depth+1),(z1,depth+1)])
  if num%2000==0:print('progress',rb,num,len(stack),accepted,pruned,round(time.time()-start,1),flush=True)
 return {'certified':True,'rho':rb,'nodes':num,'accepted_boxes':accepted,'pruned_boxes':pruned,'maxdepth':depthmax,'lower_bounds':{k:v/SCALE for k,v in lower.items()},'lower_bounds_scaled_integers':lower,'scale':SCALE,'projection_halfwidths':[h/SCALE for h in half],'density_interval':rho.pair(),'pi_interval':PI.pair(),'seconds':time.time()-start}

if __name__=='__main__':
 rb=sys.argv[1] if len(sys.argv)>1 else '0.1';print('params',m,rho,PV,PG,PC,flush=True)
 result=verify(rb);print(json.dumps(result,indent=2),flush=True)
 from pathlib import Path
 out=Path(sys.argv[2]) if len(sys.argv)>2 else Path(__file__).with_name('certificate_'+rb+'.json')
 out.write_text(json.dumps(result,indent=2))
