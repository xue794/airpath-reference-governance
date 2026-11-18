"""Independent evaluation of the frozen confirmation matrix on saved states."""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.spatial.transform import Rotation

H = Path(__file__).resolve().parents[1]
ROOT = H.parent
BASE = ROOT / 'decel_d0_20261108'
OUT = H / 'results'
sys.path.insert(0, str(BASE / 'analysis'))
from evaluate import tube, span, first, kkt_trace

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def mean_on(trace, col, lo, hi):
    part = span(trace, lo, hi)
    if len(part) < 2 or abs(part.time.iloc[0]-lo)>1e-9 or abs(part.time.iloc[-1]-hi)>1e-9:
        return np.nan
    return float(np.trapezoid(part[col], part.time)/(hi-lo))

def wind(t, gain):
    q=np.clip(t/3,0,1)
    s=q**3*(10-15*q+6*q*q)
    sd=10*q*q*(1-q)**2
    w0=np.column_stack([6+2*np.sin(.25*t),-3+1.5*np.sin(.4*t+.5),.8*np.sin(.2*t)])
    wd0=np.column_stack([.5*np.cos(.25*t),.6*np.cos(.4*t+.5),.16*np.cos(.2*t)])
    return gain*s[:,None]*w0, gain*(sd[:,None]*w0+s[:,None]*wd0)

def velocity(speed, gamma, chi):
    return (speed*np.array([np.cos(gamma)*np.cos(chi),np.cos(gamma)*np.sin(chi),-np.sin(gamma)])).T

def check_states(d, s, case, cfg):
    t=s['Ts']; x=s['Xs']; r=s['C']['R']; o=s['O']; prm=cfg['cases'][case]; xf=x[:,19:25]
    assert (o['Course'],o['SpeedRise'],o['SpeedStart'],o['Vlow'],o['Horizon'],r['rhoB'])==(0,10,3,prm['target'],35,.1)
    assert np.max(abs(d.time-t))<1e-12 and np.all(np.diff(t)>0)
    rot=Rotation.from_quat(x[:,[7,8,9,6]])
    ground=rot.apply(x[:,3:6]); w,wd=wind(t,prm['wind_gain']); air=ground-w
    gv=np.linalg.norm(ground[:,:2],axis=1); V=np.linalg.norm(air,axis=1)
    ga=np.arctan2(-air[:,2],np.linalg.norm(air[:,:2],axis=1));chi=np.arctan2(air[:,1],air[:,0])
    eg=ga-xf[:,1]; ec=np.arctan2(np.sin(chi-xf[:,2]),np.cos(chi-xf[:,2])); ev=V-xf[:,0]
    ground_ref=velocity(*xf[:,:3].T)+x[:,16:19]
    q=np.clip((t-3)/10,0,1); vreq=r['Vb']+(prm['target']-r['Vb'])*q**3*(10-15*q+6*q*q)
    for a,b in [(d.vground,gv),(d.vref_ground,np.linalg.norm(ground_ref[:,:2],axis=1)),(d.vrequest,vreq),(d.airSpeed,V),(d.refAirSpeed,xf[:,0]),(d.eV,ev),(d.egamma_deg,np.rad2deg(eg)),(d.echi_deg,np.rad2deg(ec)),(d[['windN','windE','windD']].to_numpy(),w),(d[['windDotN','windDotE','windDotD']].to_numpy(),wd),(d[['whatN','whatE','whatD']].to_numpy(),x[:,16:19])]:
        assert np.max(abs(a-b))<1e-9
    m=r['P']['m']; g=r['g']; th=r['P']['theta']; gains=s['C']['A']
    X=m*(xf[:,3]-gains['kV']*ev+g*np.sin(ga))
    ag=V*(xf[:,4]-gains['kGamma']*eg)+g*np.cos(ga)
    ac=V*np.cos(ga)*(xf[:,5]-gains['kChi']*ec)
    N=m*np.hypot(ag,ac)
    def M(alpha):
        K=.5*r['rho']*r['P']['S']*V**2
        cd=np.polynomial.polynomial.polyval(alpha,th[[0,1,4,7,9]])
        cl=np.polynomial.polynomial.polyval(alpha,th[[15,16,20,21,22]])
        return K*cl+(K*cd+X)*np.tan(alpha)
    assert max(np.max(abs(N-M(r['alphaLo'])-d.fbLo)),np.max(abs(M(r['alphaHi'])-N-d.fbHi)))<1e-9
    assert np.max(abs(d[['refLo','refHi','fbLo','fbHi']].to_numpy()-s['D'][:,12:16]))<1e-8
    P=[[Fraction(float(v)) for v in row] for row in r['PB']]
    h=np.array([tube(float(ti),xi,r,P) for ti,xi in zip(t,xf)])
    assert np.max(abs(h-d.hmove))<1e-10
    d['h_exact']=h; d['tube']=h>=0
    d['box']=(abs(ev)<=.125)&(abs(eg)<=np.deg2rad(.45))&(abs(ec)<=np.deg2rad(.35))
    d['EG']=d.tube & d.box; d['J']=d.EG & (d.theta_sum_deg<=.70)
    eligible=d.time.between(14,20)&d.EG
    assert (d[['normalMargin','bankMargin','fbHi','fbLo']][eligible]>0).all().all()
    assert (d.refHi[eligible]-d.fbHi[eligible]<10).all() and (d.refLo[eligible]-d.fbLo[eligible]<10).all()
    assert d.Traw[d.time.between(14,20)&d.J].between(0,60,inclusive='neither').all()
    d['err']=abs(gv-vreq);d['wind_norm']=np.linalg.norm(w,axis=1)
    d['wind_est_error']=np.linalg.norm(w-x[:,16:19],axis=1);d['old_ratio']=s['D'][:,16]

def event(folder,s,d):
    st=s['Summary'];last=d.iloc[-1]; status=st['status']
    if status=='completed':
        assert last.time==35
        return status,35,None
    if status=='interface_error_at_solver_trial':
        assert st['error_id']=='rg_raw_command_chart:singular',st
        ch=json.loads((folder/'chart_event.json').read_bytes());a,b=ch['runs']
        assert a['index']==b['index'] and abs(a['time_s']-b['time_s'])<1e-8
        cause=['speed_chart','horizontal_speed_chart'][b['index']-1]; ts=b['time_s'];g=b['gaps']
        data=dict(residual=b['Vchart']-s['C']['R']['Vc_min'] if b['index']==1 else b['Vhchart']-s['C']['R']['Vhc_min'],
            chartV=b['Vchart'],chartVh=b['Vhchart'],Tref=b['Tref'],Traw=b['Traw'],Tact=b['Tact'],last_saved=last.time,loc_drift_s=abs(a['time_s']-ts))
    else:
        assert status=='terminal_original_root_boundary',st
        idx=int(st['event_indices'])-1; names=['refLo','refHi','fbLo','fbHi']
        cause=['ref_lower_root','ref_upper_root','feedback_lower_root','feedback_upper_root'][idx]
        g=last[names].to_numpy(float);ts=last.time
        data=dict(residual=g[idx],chartV=last.chartV,chartVh=last.chartVh,Tref=last.Tref,Traw=last.Traw,Tact=last.Tact,last_saved=last.time)
    assert abs(data['residual'])<1e-7 and min(g)>-1e-7
    return cause,ts,dict(**data,**dict(zip(['refLo','refHi','fbLo','fbHi'],g)),cause=cause,event_s=ts)

def assess(sc,arm,step,s,d,cause,et,cover,health):
    kw=span(d,14,20); late=span(d,25,35)
    complete=cause=='completed'
    saving=float(d.time.iloc[-1])
    row=dict(case=sc,arm=arm,step=step,cause=cause,end_time=et,last_saved=saving,n=len(d),first_saved=d.time.iloc[0],
        MAE_3_16=mean_on(d,'err',3,16),position_RMS35=np.sqrt(np.mean(d.positionError**2)) if complete else np.nan,
        peak_pos_until_stop=d.positionError.max(),peak_position_s=d.time.iloc[d.positionError.argmax()],
        min_fb_hi_window=kw.fbHi.min() if saving>=20 else np.nan,min_fb_lo_window=kw.fbLo.min() if saving>=20 else np.nan)
    for k in ['vground','vref_ground','airSpeed','refAirSpeed']:
        row[k+'_25_35']=mean_on(d,k,25,35)
    row['plateau_bias']=row['vground_25_35']-s['O']['Vlow']
    for k in ['Tref','Traw','Tact','alphaActual_deg','refHi','refLo','fbHi','fbLo']:
        row[k+'_min']=d[k].min();row[k+'_max']=d[k].max()
    row['fbHi_min_time']=float(d.time.iloc[d.fbHi.argmin()])
    inter=d.change>1e-7
    row.update(interventions=inter.sum(),backup_support_on_interventions=(inter&(d.support18!=0)).sum(),endpoint_support=(d.support14.astype(bool)|d.support15.astype(bool)).sum(),
        old_out_first=first(d,d.old_ratio>1),old_out_samples=(d.old_ratio>1).sum(),
        low_out_first=first(d,(d.old_ratio>1)&(d.refAirSpeed<s['C']['R']['Vb']-.606)),
        low_plateau_samples=((late.old_ratio>1)&(late.refAirSpeed<s['C']['R']['Vb']-.606)).sum(),plateau_samples=len(late),
        KKT_max=kkt_trace(d),wind_norm_max=d.wind_norm.max(),wind_est_error_max=d.wind_est_error.max())
    assert row['KKT_max']<1e-8
    for label,y in [('saved',d),('critical',kw)]:
        if y.empty:continue
        cover.append(dict(case=sc,arm=arm,step=step,window=label,start=y.time.iloc[0],end=y.time.iloc[-1],n=len(y),tube=y.tube.sum(),airBox=y.box.sum(),EG=y.EG.sum(),J=y.J.sum(),
            tube_h_min=y.h_exact.min(),max_abs_eV=y.eV.abs().max(),max_abs_eg_deg=y.egamma_deg.abs().max(),max_abs_ec_deg=y.echi_deg.abs().max(),theta_sum_max_deg=y.theta_sum_deg.max(),
            min_kappa=y.kappa.min(),first_out_box=first(y,~y.box),first_out_tube=first(y,~y.tube),first_out_J=first(y,~y.J),
            lower_loss_max_N=(y.refLo-y.fbLo).max(),upper_loss_max_N=(y.refHi-y.fbHi).max()))
    ar=dict(case=sc,arm=arm,step=step,chartV_min=d.chartV.min(),chartVh_min=d.chartVh.min())
    for k in ['surfaceAmp','surfaceRate','thrustAmp','thrustRate','surfaceCmdAtBound','thrustCmdClip']:
        ar[k+'_samples']=int(d[k].sum());ar[k+'_approx_s']=np.trapezoid(d[k],d.time)
    health.append(ar)
    return row

def pair(case,label,a,b,key,sign,limit,floor,strict=False):
    x=a[key];y=b[key];effect=sign*(x-y)
    gate=5*max(abs(x['010']-x['005']),abs(y['010']-y['005']),floor)
    finite=np.isfinite(effect).all() and np.isfinite(gate)
    resolved=finite and np.sign(effect['005'])==np.sign(effect['010']) and effect.abs().min()>gate
    scale=resolved and (effect.min()>limit if strict else effect.min()>=limit)
    return dict(case=case,name=label,coarse=effect['010'],fine=effect['005'],drift_gate=gate,resolved=bool(resolved),meets_scale=bool(scale),research_scale=limit)

def main():
    cfg=json.loads((H/'protocol.json').read_bytes())
    for f,sha in cfg['inputs'].items(): assert digest(ROOT/f)==sha,f
    # The baseline package also pins the D0 outputs used for cross-case comparisons.
    base_manifest=json.loads((BASE/'release_manifest.json').read_bytes())
    for f in ['results/run_metrics.csv','results/decision.json','analysis/evaluate.py','results/contrasts.csv']:
        assert digest(BASE/f)==base_manifest[f]
    reg=ROOT/'decel_region_20261108'
    assert digest(reg/'results/decision.json')==json.loads((reg/'release_manifest.json').read_bytes())['results/decision.json']
    rows=[];cover=[];health=[];events=[];matched=[];free=[];common=[];traces={}; states={}
    for sc,arms in cfg['arms'].items():
        for step in ['010','005']:
            for arm in arms:
                folder=H/'runs'/f'{sc}_{arm}_{step}'
                s=loadmat(folder/(arm+'.mat'),simplify_cells=True)
                assert s['C']['MaxStep']==int(step)*.001
                d=pd.read_csv(folder/'replay.csv',float_precision='round_trip')
                check_states(d,s,sc,cfg);cause,et,e=event(folder,s,d)
                if e:events.append(dict(case=sc,arm=arm,step=step,**e))
                rows.append(assess(sc,arm,step,s,d,cause,et,cover,health))
                traces[sc,arm,step]=d;states[sc,arm,step]=s['Xs']
                if arm=='N10':
                    z=pd.read_csv(folder/'free_replay.csv')
                    assert len(z)==len(d) and np.max(abs(d.time-z.time))<1e-10 and z[['T_KKT','Am_KKT']].max().max()<1e-8
                    free.append(dict(case=sc,step=step,n=len(z),changed_T=(z.NT_delta>1e-7).sum(),changed_Am=(z.NA_delta>1e-7).sum(),first_T=first(z,z.NT_delta>1e-7),first_Am=first(z,z.NA_delta>1e-7),max_delta_T=z.NT_delta.max(),max_delta_Am=z.NA_delta.max(),violated_hi_by_T=(z.removed_phi15 < -1e-7).sum(),min_deleted_phi_hi=z.removed_phi15.min()))
            f=traces[sc,'F',step]
            for arm in ['D','A']:
                b=traces[sc,arm,step]
                i,j=np.intersect1d(f.time,b.time,return_indices=True)[1:]
                delta=f[arm+'_scaledDelta']; kkt=f['F'+arm+'KKT'].max(); assert kkt<1e-8
                matched.append(dict(case=sc,step=step,substitute=arm,n=len(f),changed=(delta>1e-7).sum(),max_scaled_delta=delta.max(),first_change=first(f,delta>1e-7),
                    common_states=len(i),same_grid=np.array_equal(f.time,b.time),max_abs_36state_diff=np.max(abs(states[sc,'F',step][i]-states[sc,arm,step][j])),
                    replay_KKT=kkt,min_deleted_phi14_by_D=f.deletedPhi14.min(),min_deleted_phi15_by_D=f.deletedPhi15.min()))
            pairs=[('F','O'),('F','D'),('F','A')]
            if sc=='C1':pairs.extend([('N10','T'),('Am','T'),('N10','Am')])
            for aa,bb in pairs:
                a=traces[sc,aa,step]; b=traces[sc,bb,step]
                y=a.merge(b,on='time',suffixes=('_a','_b'),validate='one_to_one')
                w=span(y,14,20)
                common.append(dict(case=sc,step=step,arm_a=aa,arm_b=bb,n=len(y),common_from=y.time.iloc[0],common_to=y.time.iloc[-1],
                    pos_rms_a=np.sqrt(np.mean(y.positionError_a**2)),pos_rms_b=np.sqrt(np.mean(y.positionError_b**2)),critical_from=w.time.min(),critical_to=w.time.max(),
                    min_upper_a=w.fbHi_a.min(),min_upper_b=w.fbHi_b.min(),last_upper_a=w.fbHi_a.iloc[-1] if len(w) else np.nan,last_upper_b=w.fbHi_b.iloc[-1] if len(w) else np.nan))
    metrics=pd.DataFrame(rows); contrasts=[]
    for sc in cfg['cases']:
        def arms(arm):return metrics[(metrics.arm==arm)&(metrics['case']==sc)].set_index('step')
        contrasts.append(pair(sc,'O_minus_F_speed',arms('F'),arms('O'),'MAE_3_16',-1,.606,1e-4,True))
        for other in ['D','A']:
            contrasts.append(pair(sc,'F_minus_'+other+'_upper',arms('F'),arms(other),'min_fb_hi_window',1,1,1e-7))
        if sc=='C1':contrasts.append(pair(sc,'N10_minus_Am_upper',arms('N10'),arms('Am'),'min_fb_hi_window',1,1,1e-7))
    base=pd.read_csv(BASE/'results/run_metrics.csv',dtype={'step':str}).query("arm=='F'").set_index('step').rename(columns={'governed_25_35':'vref_ground_25_35'})
    fr={sc:metrics[(metrics['case']==sc)&(metrics.arm=='F')].set_index('step') for sc in cfg['cases']}; fr['D0']=base
    diff=[]
    for a,b in [('D0','C1'),('C1','C2'),('C2','C3')]:
        for key in ['MAE_3_16','vground_25_35','vref_ground_25_35','min_fb_hi_window','position_RMS35']:
            # 5x arm drift floor applies to speed/force; position only reported as descriptive.
            x=fr[b][key];y=fr[a][key];z=x-y
            floor=1e-7 if key=='min_fb_hi_window' else 1e-4
            gate=5*max(abs(x['010']-x['005']),abs(y['010']-y['005']),floor) if key!='position_RMS35' else np.nan
            diff.append(dict(from_case=a,to_case=b,metric=key,coarse=z['010'],fine=z['005'],drift_gate=gate,resolved=bool(z.prod()>0 and z.abs().min()>gate) if np.isfinite(gate) else None))
    for name,data in [('run_metrics',metrics),('coverage',cover),('actuators',health),('events',events),('matched_replay',matched),('free_replay',free),('common_windows',common),('contrasts',contrasts),('cross_case_F',diff)]:
        pd.DataFrame(data).to_csv(OUT/(name+'.csv'),index=False,encoding='utf-8-sig')
    return metrics,traces
if __name__=='__main__':
    main()
