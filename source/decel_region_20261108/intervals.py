"""Vector intervals for binary64 inputs. Every elementary operation is widened."""
import numpy as np
from math import factorial
def down(x):return np.nextafter(x,-np.inf)
def up(x):return np.nextafter(x,np.inf)
class Box:
    __array_priority__=1000
    def __init__(self,l,h=None):
        self.l=np.asarray(l,dtype=float);self.h=self.l if h is None else np.asarray(h,dtype=float)
    def __add__(self,b):
        b=box(b);return Box(down(self.l+b.l),up(self.h+b.h))
    __radd__=__add__
    def __neg__(self):return Box(-self.h,-self.l)
    def __sub__(self,b):return self+-box(b)
    def __rsub__(self,b):return box(b)+-self
    def __mul__(self,b):
        b=box(b);a,c,d,e=self.l*b.l,self.l*b.h,self.h*b.l,self.h*b.h
        return Box(down(np.minimum(np.minimum(a,c),np.minimum(d,e))),up(np.maximum(np.maximum(a,c),np.maximum(d,e))))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=box(b);assert np.all(b.l>0)
        return self*Box(down(1/b.h),up(1/b.l))
    def sq(self):
        a=self.l*self.l;b=self.h*self.h
        return Box(np.maximum(0,down(np.where((self.l<=0)&(self.h>=0),0,np.minimum(a,b)))),up(np.maximum(a,b)))
    def sqrt(self):
        assert np.all(self.l>=0)
        return Box(np.maximum(0,down(np.sqrt(self.l))),up(np.sqrt(self.h)))
    def cap(self,b):return Box(np.maximum(self.l,b.l),np.minimum(self.h,b.h))
def box(a):return a if isinstance(a,Box) else Box(a)
def poly(c,x):
    y=box(c[-1])
    for a in c[-2::-1]:y=y*x+a
    return y
def point_trig(x,cosine):
    # For |x| <= 0.7 the first omitted Taylor term is < 0.7^16/16! < 2e-16.
    x=box(x);assert np.max(np.maximum(abs(x.l),abs(x.h)))<=.7
    c=[box((-1.)**j)/float(factorial(2*j+(not cosine))) for j in range(8)]
    y=poly(c,x.sq());y=y if cosine else y*x
    return y+Box(-2e-16,2e-16)
def trig(x):
    # Taylor on an interval (no endpoint monotonicity assumptions).
    return point_trig(x,False),point_trig(x,True)
