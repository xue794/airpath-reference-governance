"""Reconstruct every terminal leaf without using the producer's traversal.

Independently check complete dyadic partitions, exact pruning, and positive
lower bounds. Re-evaluate with monotone endpoint Taylor trigonometry (including
crossing zero for cosine), in addition to the interval polynomial producer.
"""
import csv,gzip,json,sys,time
import numpy as np
from model import OUT,geometry,save,lock
import kernel,centered,joint,intervals
from intervals import Box,box,down,up
from math import factorial

def alternate_trig(x):
    x=box(x)
    def point(z):
        z=box(z);z2=z.sq()
        def horner(cos):
            a=box(0)
            for k in range(9,-1,-1):
                den=int(factorial(2*k+(not cos)))
                inv=box((-1.)**k)/Box(down(float(den)),up(float(den)))
                a=inv-a*0+a*z2
            return (a if cos else a*z)+Box(-1e-19,1e-19)
        return horner(False),horner(True)
    # |x|<=0.7; degree 18/19 Taylor remainders < 0.7^20/20! < 1e-19.
    assert np.all(x.l>=-.7) and np.all(x.h<=.7)
    s0,c0=point(x.l);s1,c1=point(x.h)
    return Box(s0.l,s1.h),Box(np.minimum(c0.l,c1.l),np.where((x.l<=0)&(x.h>=0),1,np.maximum(c0.h,c1.h)))

def verify(mode):
    key='22_0.01_'+mode
    report=json.loads((OUT/(key+'.json')).read_bytes());assert report['certified']
    results={}
    for impl in ['original','alternative']:
        tr=intervals.trig if impl=='original' else alternate_trig
        kernel.trig=centered.trig=joint.trig=tr
        ss=report['slabs'];assert ss[0][1]==(0 if mode.startswith('nominal') else 14)
        assert ss[-1][2]==(55 if mode.startswith('nominal') else 20) and all(a[2]==b[1] for a,b in zip(ss,ss[1:]))
        proj=geometry(.01)['projection']
        nodes={}
        with gzip.open(OUT/(key+'_leaves.csv.gz'),'rt') as src:
            for row in csv.DictReader(src):
                idx=(int(row['slab']),int(row['path']));assert idx not in nodes;nodes[idx]=row['class']
        assert len(nodes)==report['nodes']
        stack=[(j,1,np.array([s[1],*(-np.array(proj))]),np.array([s[2],*proj])) for j,s in enumerate(ss)]
        batches={};count=0;outside=0;accepted=0;worsts={};per_slab={}
        def flush(phase,group):
            lower=np.stack([x[2] for x in group]);upper=np.stack([x[3] for x in group])
            z=[Box(a,b) for a,b in zip(lower.T,upper.T)]
            if mode=='nominal_focused':v=centered.sharpen(z[0],z[1:],22,phase,.01)
            elif mode=='joint':v=joint.joint(z[0],z[1:],22,phase,.01)
            else:v=kernel.feedback(z[0],z[1:],22,phase,uniform=mode.startswith('uniform'))
            for k,b in v.items():
                if mode in ('directional_focused','uniform_focused') and k.startswith('rawT'):continue
                if mode.startswith('uniform') and k.startswith('pay'):continue
                bnd=np.broadcast_to(b.l,(len(group),));j=int(np.argmin(bnd));val=float(bnd[j]);assert val>0,(mode,k,val,group[j][:2])
                if k not in worsts or val<worsts[k]['lower']:
                    worsts[k]=dict(lower=val,slab=group[j][0],path=group[j][1],time_and_offsets_low=lower[j].tolist(),time_and_offsets_high=upper[j].tolist())
                if impl=='original':
                    for i,xx in enumerate(group):
                        jsl=xx[0];per_slab.setdefault(jsl,{})[k]=min(per_slab.get(jsl,{}).get(k,np.inf),float(bnd[i]))
        while stack:
            sl,p,l,h=stack.pop();c=nodes.pop((sl,p));count+=1
            if c=='accepted':
                ph=ss[sl][0];batches.setdefault(ph,[]).append((sl,p,l,h));accepted+=1
                if len(batches[ph])>=2048:flush(ph,batches[ph]);batches[ph]=[]
            elif c=='outside':
                assert kernel.domain([Box(a,b) for a,b in zip(l[1:],h[1:])]).l>.01;outside+=1
            else:
                axis=int(c);assert 0<=axis<=5;mid=(l[axis]+h[axis])*.5
                L=l.copy();H=h.copy();L[axis]=mid;H[axis]=mid
                stack.extend([(sl,2*p,l,H),(sl,2*p+1,L,h)])
        for ph,group in batches.items():
            if group:flush(ph,group)
        assert not nodes and accepted==report['accepted'] and outside==report['pruned']
        if impl=='original':
            for k,v in worsts.items():assert abs(report['minima'][k]['lower']-v['lower'])<1e-10
        results[impl]=dict(nodes=count,accepted=accepted,pruned=outside,worsts=worsts)
        print(mode,impl,count,flush=True)
    save(OUT/(key+'_review.json'),results)
if __name__=='__main__':lock();verify(sys.argv[1])
