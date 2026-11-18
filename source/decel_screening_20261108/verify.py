"""Independent finite checks of numerical transcription and witness membership."""
from itertools import product
import json,sys
import numpy as np
import pandas as pd
from scipy.io import loadmat
from freeze import ROOT,RESULT,OLD,lock

def run():
    lock(); out={};a=loadmat(RESULT/'affine_export.mat',simplify_cells=True)
    cp=loadmat(OLD/'N10.mat',simplify_cells=True)['C'];R=cp['R'];g=R['g'];m=R['P']['m']
    assert R['DeltaN']==10 and R['rhoB']==.1 and m>0 and cp['A']['kChi']>0
    sets=dict(zip(a['labels'],a['sets']));lp=pd.read_csv(RESULT/'historical_lp.csv')
    scales=np.array([2,np.deg2rad(8),np.deg2rad(12)])
    for case,table in lp.groupby('case'):
        v=sets[case];u=table[['uV','uGamma','uChi']].values
        res=v['c']+np.einsum('nij,nj->ni',v['G'],u)
        assert np.allclose(res,table[['lp_'+n for n in v['names']]],atol=1e-8,rtol=1e-13)
        units=np.maximum(1,np.linalg.norm(v['G']*scales,axis=2))
        assert np.min(res/units)>0 and table.lp_feasible.all()
        out['LP_'+case+'_min_normalized_margin']=float(np.min(res/units))
    old=loadmat(OLD/'T.mat',simplify_cells=True)
    assert old['Summary']['event_indices']==4 and old['Ts'][-1]<16
    actual=pd.read_csv(RESULT/'actual_feedback.csv')
    emax=actual[actual.case=='N10'][['eV_mps','egamma_deg','echi_deg']].abs().max().values
    bins=np.array([.025,.05,.05]);assert np.allclose(np.ceil(emax/bins)*bins,[.125,.45,.35])
    out['empirical_box_source_max']=emax.tolist()
    out['common_prefix_ground_inputs_identical']=(np.asarray(a['prefix'])==0).all().item()
    A=cp['A'];theta=R['P']['theta']; cl=theta[[15,16,20,21,22]]; cd=theta[[0,1,4,7,9]]
    angles=[R['alphaLo'],R['alphaHi']]
    errgrid=np.array(list(product([-1,0,1],repeat=3)),dtype=float).T[:,:,None]
    primary=pd.read_csv(RESULT/'feedback_boxes.csv'); sections=pd.read_csv(RESULT/'feedback_section_probe.csv')
    checks=0;largest=0.
    for (case,box),sub in pd.concat((primary,sections),ignore_index=True).groupby(['case','box']):
        st=sets[case]
        if case.startswith('section'):ix=sub.sample_index.astype(int)
        else:
            ix=np.rint(sub.time_s*20).astype(int) if case.startswith('center') else np.arange(len(sub))
        xf=st['x'][ix].T
        assert np.allclose(st['t'][ix],sub.time_s,atol=1e-11,rtol=0)
        widths=[.125,np.deg2rad(.45),np.deg2rad(.35)] if box=='directional' else [38.1*np.deg2rad(.45),np.deg2rad(.45),np.deg2rad(.45)]
        ev,eg,ec=np.array(widths)[:,None,None]*errgrid
        v=xf[0]+ev;ga=xf[1]+eg
        X=m*(xf[3]-A['kV']*ev+g*np.sin(ga))
        N=m*np.hypot(v*(xf[4]-A['kGamma']*eg)+g*np.cos(ga),v*np.cos(ga)*(xf[5]-A['kChi']*ec))
        K=.5*R['rho']*R['P']['S']
        for j,(end,alpha) in enumerate(zip(['lower','upper'],angles)):
            lift=K*v*v*np.polynomial.polynomial.polyval(alpha,cl); drag=K*v*v*np.polynomial.polynomial.polyval(alpha,cd)
            M=lift+(X+drag)*np.tan(alpha); gap=(N-M)*(1 if j==0 else -1)
            loss=gap[13]-gap
            delta=float(max(np.max(loss.max(axis=0)-sub[end+'_loss_sup_ub_N']),
                           np.max(sub[end+'_fb_inf_lb_N']-gap.min(axis=0))))
            largest=max(largest,delta);assert delta<1e-9
            assert np.all(sub[end+'_loss_sup_lb_N']<=sub[end+'_loss_sup_ub_N']+1e-10)
            assert np.max(abs(gap[13]-sub[end+'_ref_N']))<1e-9
        checks+=len(sub)*errgrid.shape[1]
    out.update(feedback_lattice_combinations=checks,max_binary64_difference_outside_enclosure_N=largest,
               directional_rounding_pass=True,affine_dataset_count=len(sets))
    W=json.loads((RESULT/'candidate_witnesses.json').read_bytes())
    for case,wits in W.items():
        st=sets[case]
        for w in wits.values():
            offset=np.array(w['offset']);val=offset@R['PB']@offset
            assert val<st['radius'] and w['reserve_margin_N'][1]<0
            ix=np.where((np.linalg.norm(st['offset']-offset,axis=1)<1e-13)&(st['t']==w['time_s']))[0]
            assert len(ix)==1
            assert abs(st['T'][ix[0]]-np.mean(w['Tref_N']))<1e-10
    out['strictly_rejected_candidates']=list(W);out['status']='passed'
    (RESULT/'verification.json').write_text(json.dumps(out,indent=2),encoding='utf8')
    print(json.dumps(out))
if __name__=='__main__':run()
