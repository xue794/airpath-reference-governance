"""Independent six-DOF replay checks and the prespecified D0 comparisons."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
H=Path(__file__).resolve().parents[1];ROOT=H.parent;OUT=H/'results';OUT.mkdir(exist_ok=True)
ARMS=['O','F','D','A','T','N10','Am']
def put(name,d):
    (OUT/name).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
def span(d,a,b):
    return d[(d.time>=a-1e-10)&(d.time<=b+1e-10)]
def avg(d,key):
    return float(np.trapezoid(d[key],d.time)/(d.time.iloc[-1]-d.time.iloc[0]))
def first(d,mask):
    t=d.time[mask];return None if t.empty else float(t.iloc[0])
def tube(t,x,R,P):
    t=F(t);v0=F(float(R['Vb']));v1=F(125,4)
    q=min(F(1),max(F(0),(t-3)/22))
    v=v0+(v1-v0)*q**3*(10-15*q+6*q*q);n=(v1-v0)*30*q*q*(1-q)**2/22
    a,b,c,d,e=[F(float(x[k])) for k in [0,1,3,4,5]];a-=v;c-=n
    z=P[0][0]*a*a+2*P[0][2]*a*c+P[2][2]*c*c+P[1][1]*b*b+2*P[1][3]*b*d+P[3][3]*d*d+P[4][4]*e*e
    return float(F(.01)-z)
def kkt_trace(d):
    return float(d[['primal','stationarity','dualViolation','complementarity']].max().max())
def run():
    cfg=json.loads((H/'protocol.json').read_bytes())
    for path,h in cfg['inputs'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    records=[];cert=[];prefix=[];same=[];traces={};data={};exits=[];health=[]
    for step in ['010','005']:
        for arm in ARMS:
            folder=H/'runs'/(arm+'_'+step)
            s=loadmat(folder/(arm+'.mat'),simplify_cells=True)
            d=pd.read_csv(folder/'replay.csv',float_precision='round_trip');t=s['Ts'];x=s['Xs'];raw=s['D'];C=s['C'];R=C['R'];st=s['Summary']
            assert s['O']['Course']==0 and s['O']['SpeedRise']==8 and s['O']['Vlow']==31 and R['rhoB']==.1
            assert np.max(abs(d.time-t))<1e-12 and np.all(np.diff(t)>0) and kkt_trace(d)<1e-8
            vel=Rotation.from_quat(x[:,[7,8,9,6]]).apply(x[:,3:6])
            gv=np.hypot(vel[:,0],vel[:,1]);V=np.linalg.norm(vel,axis=1);ga=np.arctan2(-vel[:,2],gv);chi=np.arctan2(vel[:,1],vel[:,0])
            xf=x[:,19:25];eg=ga-xf[:,1];ec=np.arctan2(np.sin(chi-xf[:,2]),np.cos(chi-xf[:,2]));ev=V-xf[:,0]
            dt=np.clip((t-3)/8,0,1);vreq=R['Vb']+(31-R['Vb'])*dt**3*(10-15*dt+6*dt**2)
            assert max(np.max(abs(d.vground-gv)),np.max(abs(d.vrequest-vreq)),np.max(abs(d.eV-ev)),np.max(abs(d.egamma_deg-np.rad2deg(eg))),np.max(abs(d.echi_deg-np.rad2deg(ec))))<1e-9
            # Independent force evaluation including the unrestricted terminal event residual.
            m=R['P']['m'];g=R['g'];theta=R['P']['theta'];A=C['A']
            X=m*(xf[:,3]-A['kV']*ev+g*np.sin(ga));ag=V*(xf[:,4]-A['kGamma']*eg)+g*np.cos(ga);ac=V*np.cos(ga)*(xf[:,5]-A['kChi']*ec)
            N=m*np.hypot(ag,ac)
            def M(alpha):
                K=.5*R['rho']*R['P']['S']*V**2
                cd=np.polynomial.polynomial.polyval(alpha,theta[[0,1,4,7,9]])
                cl=np.polynomial.polynomial.polyval(alpha,theta[[15,16,20,21,22]])
                return K*cl+(K*cd+X)*np.tan(alpha)
            assert np.max(abs(N-M(0)-d.fbLo))<1e-9 and np.max(abs(M(R['alphaHi'])-N-d.fbHi))<1e-9
            assert np.max(abs(d[['refLo','refHi','fbLo','fbHi']].values-raw[:,12:16]))<1e-8
            p=[[F(float(v)) for v in row] for row in R['PB']]
            h=np.array([tube(float(ti),xi,R,p) for ti,xi in zip(t,xf)])
            assert np.max(abs(h-d.hmove))<1e-10
            d['exact_move_h']=h; d['old_VB_ratio']=raw[:,16];d['err']=abs(gv-vreq)
            d['in_move']=h>=0;d['in_box']=(abs(ev)<=.125)&(abs(eg)<=np.deg2rad(.45))&(abs(ec)<=np.deg2rad(.35))
            d['in_EG']=d.in_move&d.in_box;d['in_J']=d.in_EG&(d.theta_sum_deg<=.70)
            d['loss_upper']=d.refHi-d.fbHi;d['loss_lower']=d.refLo-d.fbLo
            assert np.all(d[['normalMargin','bankMargin','fbHi','fbLo']][d.in_EG&d.time.between(14,20)]>0)
            assert np.all(d.Traw[d.in_J&d.time.between(14,20)].between(0,60,inclusive='neither'))
            n=len(t);traces[arm,step]=d;data[arm,step]=s;kw=span(d,14,20);win=span(d,3,16);late=span(d,25,35)
            inter=d.change>cfg['intervention_scaled_tol'];end=float(t[-1])
            et=end;cause=st['status']
            if cause=='interface_error_at_solver_trial':
                assert st['error_id']=='rg_raw_command_chart:singular'
                chart=json.loads((folder/'chart_event.json').read_bytes());a,b=chart['runs'];assert a['index']==b['index']==1 and abs(a['time_s']-b['time_s'])<1e-8
                et=b['time_s'];cause='raw_command_speed_chart'
                exits.append(dict(arm=arm,step=step,last_saved=end,t=et,event=cause,chartV=b['Vchart'],chartVh=b['Vhchart'],four_gaps=b['gaps'],Tref=b['Tref'],Traw=b['Traw'],Tact=b['Tact']))
            elif cause=='terminal_original_root_boundary':
                assert st['event_indices']==4 and abs(d.fbHi.iloc[-1])<1e-7
                cause='feedback_upper_root';exits.append(dict(arm=arm,step=step,last_saved=end,t=et,event=cause,chartV=d.chartV.iloc[-1],chartVh=d.chartVh.iloc[-1],four_gaps=d[['refLo','refHi','fbLo','fbHi']].iloc[-1].tolist(),Tref=d.Tref.iloc[-1],Traw=d.Traw.iloc[-1],Tact=d.Tact.iloc[-1]))
            else:assert cause=='completed' and end==35
            records.append(dict(arm=arm,step=step,n=n,saved_from=t[0],end_time=et,cause=cause,MAE_3_16=avg(win,'err') if end>=16 else None,
                vground_25_35=avg(late,'vground') if end==35 else None,plateau_bias=avg(late,'vground')-31 if end==35 else None,
                governed_25_35=avg(late,'vref_ground') if end==35 else None,
                position_RMS35=float(np.sqrt(np.mean(d.positionError**2))) if end==35 else None,peak_position_to_stop=d.positionError.max(),
                min_fb_hi_window=kw.fbHi.min() if end>=20 else None,min_fb_lo_window=kw.fbLo.min() if end>=20 else None,
                fbHi_to_stop=d.fbHi.min(),fbLo_to_stop=d.fbLo.min(),t_min_hi=d.time.iloc[d.fbHi.argmin()],
                min_Tref=d.Tref.min(),max_Tref=d.Tref.max(),min_Traw=d.Traw.min(),max_Traw=d.Traw.max(),min_Tact=d.Tact.min(),max_Tact=d.Tact.max(),
                min_alpha=d.alphaActual_deg.min(),max_alpha=d.alphaActual_deg.max(),interventions=inter.sum(),backup_support_interventions=(inter&(d.support18==1)).sum(),
                lower_endpoint_support=d.support14.sum(),upper_endpoint_support=d.support15.sum(),
                first_old_out=first(d,d.old_VB_ratio>1),old_out_samples=(d.old_VB_ratio>1).sum(),reg_min_h=h.min()))
            for part,sy in [('saved',d),('critical',kw)]:
                cert.append(dict(arm=arm,step=step,window=part,n=len(sy),start=sy.time.iloc[0],end=sy.time.iloc[-1],tube=sy.in_move.sum(),airBox=sy.in_box.sum(),EG_eligible=sy.in_EG.sum(),J_eligible=sy.in_J.sum(),theta_max_deg=sy.theta_sum_deg.max(),min_kappa=sy.kappa.min(),max_abs_eV=sy.eV.abs().max(),max_abs_eg_deg=sy.egamma_deg.abs().max(),max_abs_echi_deg=sy.echi_deg.abs().max(),first_outside_tube=first(sy,~sy.in_move),first_outside_box=first(sy,~sy.in_box),lower_loss_max_N=sy.loss_lower.max(),upper_loss_max_N=sy.loss_upper.max()))
            hp=dict(arm=arm,step=step,KKTmax=kkt_trace(d),chartVmin=d.chartV.min(),chartVhmin=d.chartVh.min(),surface_actual_min_rad=x[:,32:35].min(axis=0).tolist(),surface_actual_max_rad=x[:,32:35].max(axis=0).tolist())
            for k in ['surfaceAmp','surfaceRate','thrustAmp','thrustRate','surfaceCmdAtBound','thrustCmdClip']:
                hp[k+'_samples']=int(d[k].sum());hp[k+'_approx_s']=float(np.trapezoid(d[k],t))
            health.append(hp)
            if arm=='F':
                same.append(dict(step=step,n=n,changed_D=int((d.D_scaledDelta>1e-7).sum()),changed_A=int((d.A_scaledDelta>1e-7).sum()),maxdiff_D=d.D_scaledDelta.max(),maxdiff_A=d.A_scaledDelta.max(),min_deleted_phi14=d.deletedPhi14.min(),min_deleted_phi15=d.deletedPhi15.min(),max_replay_KKT=d[['FDKKT','FAKKT']].max().max()))
            if arm in ['T','N10','Am']:
                folder=(ROOT/'contribution_probe_20261106/full6dof'/('V31_C60_fine' if step=='005' else 'V31_C60')) if arm=='Am' else ROOT/'endpoint_probe_20261106/full6dof'/('V31_C60_warm_fine' if step=='005' else 'V31_C60_warm')
                old=loadmat(folder/(arm+'.mat'),simplify_cells=True);mask=t<=16;other=old['Ts']<=16
                assert np.array_equal(t[mask],old['Ts'][other]);delta=float(np.max(abs(x[mask]-old['Xs'][other])));assert delta<cfg['prefix_tol_state_abs']
                prefix.append(dict(arm=arm,step=step,n=mask.sum(),end=t[mask][-1],max_abs_36state_difference=delta))
    sums=pd.DataFrame(records)
    for name,rows in [('run_metrics',records),('coverage',cert),('matched_replay',same),('prefix',prefix),('events',exits),('actuators',health)]:
        pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False,encoding='utf-8-sig')
    contrasts=[]
    for label,a,b,metric,sign,threshold,floor in [('O_F_speed','F','O','MAE_3_16',-1,.606,.0001),('F_D_margin','F','D','min_fb_hi_window',1,1,1e-7),('F_A_margin','F','A','min_fb_hi_window',1,1,1e-7),('N10_Am_margin','N10','Am','min_fb_hi_window',1,1,1e-7)]:
        pairs=[sums[sums.arm==m].set_index('step')[metric] for m in [a,b]]
        x,y=pairs;z=sign*(x-y); gate=5*max(abs(x['010']-x['005']),abs(y['010']-y['005']),floor)
        contrasts.append(dict(name=label,coarse=z['010'],fine=z['005'],drift_gate=gate,resolved=bool(z.prod()>0 and z.abs().min()>gate),effect_scale=threshold,meets_scale=bool(z.min()>gate and z.min()>=threshold)))
    pd.DataFrame(contrasts).to_csv(OUT/'contrasts.csv',index=False,encoding='utf-8-sig')
    # Compare only matched retained grid points at early termination.
    common=[]
    for step in ['010','005']:
        for a,b in [('F','O'),('N10','T'),('Am','T'),('N10','Am'),('F','D'),('F','A')]:
            x=traces[a,step].copy();y=traces[b,step].copy()
            z=pd.merge(x,y,on='time',suffixes=('_a','_b'),validate='one_to_one')
            w=span(z,14,20)
            common.append(dict(step=step,arm_a=a,arm_b=b,n=len(z),from_s=z.time.iloc[0],to_s=z.time.iloc[-1],pos_rms_a=np.sqrt(np.mean(z.positionError_a**2)),pos_rms_b=np.sqrt(np.mean(z.positionError_b**2)),critical_start=w.time.iloc[0],critical_end=w.time.iloc[-1],min_upper_a=w.fbHi_a.min(),min_upper_b=w.fbHi_b.min()))
    pd.DataFrame(common).to_csv(OUT/'common_windows.csv',index=False,encoding='utf-8-sig')
    for (arm,step),d in traces.items():
        d.to_csv(OUT/f'{arm}_{step}_trace.csv',index=False,encoding='utf-8-sig')
    return sums,contrasts,traces
if __name__=='__main__':run()
