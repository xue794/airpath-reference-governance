"""Replay the archived backup, solve point feasibility and check nonzero sections."""
import json
import numpy as np
from scipy.io import loadmat
from scipy.optimize import linprog
from freeze import RESULT,OLD,AUDIT,lock
from screen_feedback import wp,I,pt,lo,hi,guide,ForceCheck,csvout

def dump(name,v):
    (RESULT/name).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')

def center_witness(engine,R,s,idx,d,kind):
    # Exact center at prescribed time, stored float offsets (= admissible points,
    # not the rounded MATLAB sums). Recompute e'Pe and the actual state margins.
    t=float(s['t'][idx]); old=guide(t,d);err=list(map(pt,s['offset'][idx]))
    x=[old[0]+err[0],err[1],pt(0),old[3]+err[2],err[3],err[4]]
    val=sum(err[j]*pt(R['PB'][j,k])*err[k] for j in range(5) for k in range(5))
    raw,q,X=engine.evaluate(x,(pt(0),)*3)
    slope=engine.slopes(q,X)
    assert hi(val)<float(s['radius']) and lo(raw[0])>0 and lo(raw[1])>0 and slope>0
    m=engine.model;N=m.m*I.sqrt((x[0]*x[4]+m.g*I.cos(x[1]))**2+(x[0]*I.cos(x[1])*x[5])**2)
    a0,a1=m.alo,m.ahi
    for _ in range(47):
        mid=(a0+a1)/2; a=pt(mid)
        residual=q*(wp.poly(m.cl,a)+wp.poly(m.cd,a)*I.tan(a))+X*I.tan(a)-N
        if lo(residual)>0:a1=mid
        elif hi(residual)<0:a0=mid
        else:break
    alpha=I.mpf([a0,a1]);T=(q*wp.poly(m.cd,alpha)+X)/I.cos(alpha)
    nplus=raw[1]-pt(10); margin=T-3 if kind=='T_min' else nplus
    assert hi(margin)<0 and a1-a0<1e-10
    return dict(time_s=t,test=kind,offset=list(map(float,s['offset'][idx])),VB_upper=hi(val),
                radius=float(s['radius']),reserve_margin_N=[lo(margin),hi(margin)],
                Tref_N=[lo(T),hi(T)],hNplus_N=[lo(nplus),hi(nplus)],reference_slope_lb=slope)

def run():
    lock()
    b=loadmat(RESULT/'affine_export.mat',simplify_cells=True)
    original=loadmat(AUDIT,simplify_cells=True)
    model=loadmat(OLD/'N10.mat',simplify_cells=True)
    assert np.max(abs(b['prefix']))==0 and model['C']['R']['rhoB']==.1
    engine=ForceCheck(model['C']);lp_rows=[];samples=[];mins=[];groups={};witnesses={}
    R=b['R'];sc=np.asarray(model['O']['QP']['Scale'])
    for name,s in zip(b['labels'],b['sets']):
        x=s['x'];ts=s['t'];ph=s['phi'];n=len(ts);names=list(s['names'])
        assert np.all(np.isfinite(ph)), f'root domain exit {name}'
        if name=='N10':
            assert np.array_equal(x,original['X']) and np.array_equal(ts,original['t'])
            assert np.max(abs(ph-original['Phi'][:,:,3]))<1e-10
            assert np.max(abs(s['vb']-original['VB'][:,3]))<1e-10
        fields=list(names); val=ph
        stat=np.column_stack([s['h'],s['psi'],s['hRate'],s['hGeom'],s['hForce'],s['hThrust']])
        state_names=['hVlo','hVhi','hGammaLo','hGammaHi','psiVlo','psiVhi','psiGlo','psiGhi','hRateV','hRateG','hRateChi','hNormal','hBankPlus','hBankMinus','hNminus','hNplus','hTmin','hTmax']
        for typ,mat,keys in [('row',ph,names),('state',stat,state_names)]:
            for j,key in enumerate(keys):
                idx=int(np.argmin(mat[:,j])); mins.append(dict(case=name,kind=typ,item=key,min=float(mat[idx,j]),time_s=ts[idx],
                              negatives=int(sum(mat[:,j]<-1e-9)),index=idx))
        # Cover M_alpha over the full alpha interval simultaneously at all
        # sampled states. Bounds use conservative hulls with outward arithmetic.
        V=I.mpf([float(min(x[:,0])),float(max(x[:,0]))]);nv=I.mpf([float(min(x[:,3])),float(max(x[:,3]))])
        ga=I.mpf([float(min(x[:,1])),float(max(x[:,1]))])
        m=engine.model; slope=engine.slopes(m.K*V**2,m.m*(nv+m.g*I.sin(ga)))
        assert slope>0
        good=(ph.min(axis=1)>0)&(stat.min(axis=1)>0)
        groups[name]=dict(n=n,positive_all=int(sum(good)),phi_neg=int(sum(np.min(ph,axis=1)<-1e-9)),
                       slope_lb_N_per_rad=slope,Tmin=float(min(s['T'])),Tmax=float(max(s['T'])),
                       min_hNplus=float(min(s['hForce'][:,1])))
        if name.startswith('section'):
            err=s['offset'];assert np.all(np.einsum('ni,ij,nj->n',err,R['PB'],err)<s['radius'])
            if not np.all(good):
                why={}
                for field,col,kind in [('hThrust',0,'T_min'),('hForce',1,'N_plus')]:
                    idx=int(np.argmin(s[field][:,col]))
                    if s[field][idx,col]<0:
                        why[kind]=center_witness(engine,R,s,idx,int(name[7:9]),kind)
                assert why
                witnesses[name]=why
        for i,t in enumerate(ts):
            xf=x[i]; c=s['c'][i];G=s['G'][i]
            row=dict(case=name,time_s=t,sample_index=i,V=xf[0],gamma=xf[1],chi=xf[2],nuV=xf[3],nuGamma=xf[4],nuChi=xf[5],
                uV=s['u'][i,0],uGamma=s['u'][i,1],uChi=s['u'][i,2],VB_original=s['vb'][i],Tref_N=s['T'][i],Ma_root=s['slope'][i])
            row.update(dict(zip(state_names,stat[i])));row.update({'phi_'+k:v for k,v in zip(names,ph[i])})
            samples.append(row)
            if name not in ('N10','T'):continue
            weights=np.maximum(1,np.linalg.norm(G*sc,axis=1))
            A=np.column_stack((-G*sc/weights[:,None],np.ones(17)))
            result=linprog([0,0,0,-1],A_ub=A,b_ub=c/weights,bounds=[(None,None)]*4,method='highs',
                options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
            assert result.success, (name,i,result.message)
            u=result.x[:3]*sc;res=c+G@u
            assert abs(min(res/weights)-result.x[-1])<1e-8
            # A feasible primal input is enough for the existential sample test.
            lp_rows.append(dict(case=name,time_s=t,VB_over_rho=s['vb'][i]/.1,backup_negative='|'.join(np.array(names)[ph[i]<-1e-9]),
                lp_margin_scaled=result.x[-1],lp_feasible=int(min(res/weights)>=-1e-9),uV=u[0],uGamma=u[1],uChi=u[2],
                **{'backup_'+k:v for k,v in zip(names,ph[i])},**{'lp_'+k:v for k,v in zip(names,res)}))
    csvout('recovery_rows.csv',samples);csvout('row_minima.csv',mins);csvout('historical_lp.csv',lp_rows)
    dump('candidate_witnesses.json',witnesses);dump('recovery_summary.json',groups)
if __name__=='__main__':run()
