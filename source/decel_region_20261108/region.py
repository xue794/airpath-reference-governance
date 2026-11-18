"""Complete bisection forest indexed by (initial time slab, binary heap path)."""
import argparse,gzip,csv,time
import numpy as np
from model import OUT,save,geometry,lock
from intervals import Box
from kernel import nominal,feedback,domain
def slabs(D,mode):
    if mode!='nominal':return [(1,i/4,(i+1)/4) for i in range(56,80)]
    seg=[0,3,3+D,35,55];result=[]
    for phase,(a,b) in enumerate(zip(seg,seg[1:])):
        if phase in (0,2):result.append((phase,a,b))
        else:result.extend((phase,i/4,(i+1)/4) for i in range(4*a,4*b))
    return result
def run(D,r,mode,refine=False):
    lock();key=f'{D}_{r}_{mode}';ss=slabs(D,mode);h=geometry(r)['projection']
    cap=1000000
    if refine:
        import json
        from centered import sharpen
        prior=OUT/(key+'.json')
        if prior.exists():cap-=json.loads(prior.read_bytes())['nodes']
        key+='_focused'
    if mode=='joint':
        import json
        from joint import joint
        cap=1000000-sum(json.loads((OUT/f'{D}_{r}_{k}.json').read_bytes())['nodes'] for k in ['directional','directional_focused'])
    start=time.process_time();nodes=0;accepted=pruned=0;minima={};bytime={};unresolved=[];globalnames=None
    kernels={'nominal':nominal,'directional':feedback,'uniform':feedback}
    with gzip.open(OUT/(key+'_leaves.csv.gz'),'wt',newline='') as fp:
        wr=csv.writer(fp);wr.writerow(['slab','path','class'])
        for phase in sorted({s[0] for s in ss}):
            batch=[(j,1,np.array([a,*(-np.array(h))]),np.array([b,*h])) for j,(ph,a,b) in enumerate(ss) if ph==phase]
            depth=0
            while batch:
                pending=[]
                for pos in range(0,len(batch),2048):
                    part=batch[pos:pos+2048];nodes+=len(part)
                    lo=np.stack([b[2] for b in part]);hi=np.stack([b[3] for b in part])
                    zz=[Box(l,u) for l,u in zip(lo.T,hi.T)];ex=domain(zz[1:]).l>r
                    for j in np.flatnonzero(ex):wr.writerow([*part[j][:2],'outside'])
                    pruned+=sum(ex);ids=np.flatnonzero(~ex)
                    if not len(ids):continue
                    z=[Box(b.l[ids],b.h[ids]) for b in zz]
                    if mode=='joint':vals=joint(z[0],z[1:],D,phase,r)
                    elif refine and mode=='nominal':vals=sharpen(z[0],z[1:],D,phase,r)
                    else:vals=kernels[mode](z[0],z[1:],D,phase,**({'uniform':True} if mode=='uniform' else {}))
                    names=list(vals);bounds=np.stack([np.broadcast_to(b.l,(len(ids),)) for b in vals.values()],axis=1)
                    assert np.isfinite(bounds).all()
                    # Uniform budget already has strict interior counterexamples; separately check its E/G/J.
                    required=[i for i,n in enumerate(names) if (mode!='uniform' or not n.startswith('pay')) and not (refine and mode!='nominal' and n.startswith('rawT'))]
                    ok=(bounds[:,required].min(axis=1)>0)
                    for i in np.flatnonzero(ok):
                        j=ids[i];sl,path=part[j][:2];wr.writerow([sl,path,'accepted']);accepted+=1
                        if sl not in bytime:bytime[sl]=bounds[i].copy()
                        else:bytime[sl]=np.minimum(bytime[sl],bounds[i])
                    globalnames=names
                    for i in np.flatnonzero(~ok):
                        j=ids[i];sl,path,l,u=part[j]
                        reserve=len(ss)*2 if refine or mode=='joint' else 10
                        if len(unresolved) or depth>=49 or nodes+len(batch)+len(pending)>=cap-reserve or time.process_time()-start>14380:
                            unresolved.append((sl,path));wr.writerow([sl,path,'unresolved']);continue
                        # Focus nominal dimensions on the blocking physical condition.
                        # Fixed schedule per condition is stored as an explicit axis at each split.
                        bad={names[k] for k in required if bounds[i,k]<=0}
                        order=[3,2,4,1,5,0] if mode in ('nominal','joint') else [5,4,1,2,3,0]
                        if mode=='nominal' and all('N' in n for n in bad):order=[4,1,5,2,3,0]
                        widths=(u-l)/np.array([.25,*np.multiply(h,2)])
                        if phase in (0,2):widths[0]=0
                        else:widths[0]*=(.8 if refine and mode=='nominal' else .18)
                        jax=order[int(np.argmax(widths[order]))];mid=(l[jax]+u[jax])/2
                        wr.writerow([sl,path,str(jax)])
                        u1=u.copy();u1[jax]=mid;l1=l.copy();l1[jax]=mid
                        pending.extend([(sl,2*path,l,u1),(sl,2*path+1,l1,u)])
                batch=pending;depth+=1
                if depth%4==0 and batch: print(key,phase,depth,nodes,len(batch),flush=True)
    for i,n in enumerate(globalnames):
        sl=min(bytime,key=lambda k:bytime[k][i]);minima[n]={'lower':float(bytime[sl][i]),'slab':ss[sl][1:]}
    summary=dict(mode=mode,D=D,r=r,certified=not unresolved,nodes=int(nodes),accepted=int(accepted),pruned=int(pruned),
                 unresolved=len(unresolved),cpu_s=time.process_time()-start,minima=minima,geometry=geometry(r),slabs=ss)
    save(OUT/(key+'.json'),summary)
    with open(OUT/(key+'_slabs.csv'),'w',newline='') as f:
        w=csv.writer(f);w.writerow(['slab','phase','start','end',*globalnames])
        for j,v in sorted(bytime.items()):w.writerow([j,*ss[j],*v])
    print(key,'complete',summary['certified'],nodes,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['nominal','directional','uniform','joint']);p.add_argument('r',type=float);p.add_argument('--D',type=int,default=22);p.add_argument('--focused',action='store_true');args=p.parse_args();run(args.D,args.r,args.mode,args.focused)
