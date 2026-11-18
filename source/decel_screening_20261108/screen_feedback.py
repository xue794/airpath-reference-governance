"""Section 2B. Enclose the error box at each prescribed reference point."""
import csv,itertools,math,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from freeze import HERE,ROOT,OLD,FB,RESULT,lock
sys.path.insert(0,str(ROOT/'revision_wp3_20261107'))
import interval_wp3 as wp

I=wp.I; pt=wp.point; lo=wp.low; hi=wp.high
HALF=I.mpf('0.45')*I.pi/180
BOXES={'directional':(I.mpf('0.125'),HALF,I.mpf('0.35')*I.pi/180),
       'uniform':(I.mpf('38.1')*HALF,HALF,HALF)}
SIGNS=list(itertools.product([-1,1],repeat=3))

def csvout(name,rows):
    with (RESULT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def guide(t,d):
    # Protocol uses analytic center; interval inputs enclose exact rational sample times.
    if t<=3:return (pt(38.1),pt(0),pt(0),pt(0),pt(0),pt(0))
    if t<=3+d:a,b,q,h=pt(38.1),pt(31.25),(I.mpf(str(t))-3)/d,d
    elif t<=35:a,b,q,h=pt(31.25),pt(31.25),pt(0),1
    else:a,b,q,h=pt(31.25),pt(38.1),(I.mpf(str(t))-35)/20,20
    v=a+(b-a)*q**3*(10-15*q+6*q**2);nv=(b-a)*30*q**2*(1-q)**2/h
    return (v,pt(0),pt(0),nv,pt(0),pt(0))

class ForceCheck:
    def __init__(self,C):
        self.model=wp.Evaluator(C)
        self.terms=[]
        for k in range(8):
            a=pt(self.model.alo)+(pt(self.model.ahi)-pt(self.model.alo))*I.mpf([k,k+1])/8
            cl,cd=self.model.cl,self.model.cd
            coeff=wp.poly([i*cl[i] for i in range(1,len(cl))],a)
            coeff+=wp.poly([i*cd[i] for i in range(1,len(cd))],a)*I.tan(a)+wp.poly(cd,a)/I.cos(a)**2
            self.terms.append((coeff,1/I.cos(a)**2))

    def evaluate(self,x,errs):
        m=self.model;ev,eg,ec=errs
        v=x[0]+ev;ga=x[1]+eg;q=m.K*v**2
        X=m.m*(x[3]-m.kv*ev+m.g*I.sin(ga))
        ag=v*(x[4]-m.kg*eg)+m.g*I.cos(ga)
        ac=v*I.cos(ga)*(x[5]-m.kc*ec)
        N=m.m*I.sqrt(ag**2+ac**2)
        return (N-q*m.coefflo-X*m.tanlo,q*m.coeffhi+X*m.tanhi-N),q,X

    def slopes(self,q,X):
        return min(lo(q*c+X*s) for c,s in self.terms)

    def minimum_faces(self,x,widths,end):
        # For given V,gamma, N increases with |nuChi-kChi*eChi|.
        # Reduce to the known extreme in course error, then use strict signs of
        # partial derivatives over the remaining rectangle to reduce other axes.
        m=self.model
        e=[I.mpf([-1,1])*w for w in widths]
        if end==1:
            if x[5]>=0:e[2]=-widths[2]
            elif x[5]<=0:e[2]=widths[2]
        else:
            center=x[5]/m.kc
            if center>=widths[2]:e[2]=widths[2]
            elif center<=-widths[2]:e[2]=-widths[2]
            elif center>=-widths[2] and center<=widths[2]:e[2]=center
        v=x[0]+e[0];ga=x[1]+e[1];co=I.cos(ga);si=I.sin(ga)
        ng=x[4]-m.kg*e[1];nc=x[5]-m.kc*e[2]
        ag=v*ng+m.g*co;ac=v*co*nc;r=I.sqrt(ag**2+ac**2)
        assert lo(r)>0 and lo(m.kc)>0
        dn_v=m.m*(ag*ng+ac*co*nc)/r
        dn_g=m.m*(ag*(-v*m.kg-m.g*si)-ac*v*si*nc)/r
        coeff,tan=(m.coefflo,m.tanlo) if end==0 else (m.coeffhi,m.tanhi)
        deriv=[dn_v-(2*m.K*v*coeff-m.m*m.kv*tan),dn_g-m.m*m.g*co*tan]
        for j,value in enumerate(deriv):
            vj=value if end==0 else -value
            if lo(vj)>0:e[j]=-widths[j]
            elif hi(vj)<0:e[j]=widths[j]
        return self.evaluate(x,e)[0][end]

    def examine(self,x,name):
        widths=BOXES[name];box=[I.mpf([-1,1])*h for h in widths]
        ref=self.evaluate(x,(pt(0),)*3)[0]
        fb,q,X=self.evaluate(x,box)
        slope=self.slopes(q,X)
        # The 8 actual vertices bound the maximum from below, even when an
        # extremum is interior. The interval over the entire box bounds it above.
        vert=[self.evaluate(x,[h*s for h,s in zip(widths,signs)])[0] for signs in SIGNS]
        row={'box':name,'slope_lower_N_per_rad':slope}
        for i,end in enumerate(['lower','upper']):
            # Both independent enclosures include the true global minimum.
            face=self.minimum_faces(x,widths,i)
            loss_ub=min(hi(ref[i]-fb[i]),hi(ref[i]-face));v=[lo(ref[i]-g[i]) for g in vert]
            j=int(np.argmax(v))
            row.update({f'{end}_ref_N':float(ref[i].mid),f'{end}_loss_sup_lb_N':v[j],
                        f'{end}_loss_sup_ub_N':loss_ub,f'{end}_fb_inf_lb_N':max(lo(fb[i]),lo(face)),
                        f'{end}_fb_inf_ub_N':min(hi(g[i]) for g in vert),
                        f'{end}_vertex_sign':''.join('+' if s==1 else '-' for s in SIGNS[j])})
        row['pay_10N']='yes' if max(row[f'{e}_loss_sup_ub_N'] for e in ('lower','upper'))<10 else (
            'counterexample' if max(row[f'{e}_loss_sup_lb_N'] for e in ('lower','upper'))>10 else 'unresolved')
        row['E']='yes' if slope>0 and min(row[f'{e}_fb_inf_lb_N'] for e in ('lower','upper'))>0 else (
            'counterexample' if min(row[f'{e}_fb_inf_ub_N'] for e in ('lower','upper'))<0 else 'unresolved')
        return row

def run():
    lock(); rows=[];actual=[]
    for name in ['N10','T']:
        d=loadmat(OLD/(name+'.mat'),simplify_cells=True); engine=ForceCheck(d['C'])
        old=list(csv.DictReader((FB/f'V31_C60_warm_fine_{name}_transfer.csv').open()))
        assert len(old)==len(d['Ts'])
        for idx,(t,xs,prior) in enumerate(zip(d['Ts'],d['Xs'],old)):
            x=xs[19:25]; vin=Rotation.from_quat(np.r_[xs[7:10],xs[6]]).apply(xs[3:6])
            angles=np.array([math.atan2(-vin[2],math.hypot(*vin[:2])),math.atan2(vin[1],vin[0])])
            e=np.r_[np.linalg.norm(vin)-x[0],angles[0]-x[1],math.remainder(angles[1]-x[2],2*math.pi)]
            ref=engine.model.numeric(x,np.zeros(3))[0];fb=engine.model.numeric(x,e)[0]
            expected=[float(prior[k]) for k in ['reference_lower_N','reference_upper_N','feedback_lower_N','feedback_upper_N']]
            assert abs(float(prior['time_s'])-t)<1e-10 and np.max(abs(np.r_[ref,fb]-expected))<1e-8
            assert np.max(abs(np.r_[ref,fb]-d['D'][idx,12:16]))<1e-8
            actual.append(dict(case=name,time_s=t,eV_mps=e[0],egamma_deg=math.degrees(e[1]),echi_deg=math.degrees(e[2]),
                ref_lower_N=ref[0],ref_upper_N=ref[1],fb_lower_N=fb[0],fb_upper_N=fb[1],
                loss_lower_N=ref[0]-fb[0],loss_upper_N=ref[1]-fb[1],Tref_N=d['D'][idx,9],Traw_N=d['D'][idx,10]))
            for box in BOXES:
                rows.append(dict(case=name,time_s=t,**engine.examine(tuple(map(pt,x)),box)))
    for d in [18,22,26]:
        for tick in range(280,401):
            t=tick/20
            for box in BOXES:
                rows.append(dict(case=f'center{d}',time_s=t,**engine.examine(guide(t,d),box)))
    csvout('feedback_boxes.csv',rows);csvout('actual_feedback.csv',actual)

if __name__=='__main__': run()
