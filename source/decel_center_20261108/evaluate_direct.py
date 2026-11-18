"""Check the archived and new trajectories from their physical states; evaluate the preregistered contrasts."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.spatial.transform import Rotation

H=Path(__file__).resolve().parent
PREV=H.parent/'decel_confirm_20261108'
OUT=H/'results'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ramp(q):return q**3*(10+q*(-15+6*q))
def averaged(d,col,a,b):
    q=d[d.time.between(a-1e-9,b+1e-9)]
    if len(q)<2 or not np.allclose(q.time.iloc[[0,-1]],[a,b],atol=1e-9,rtol=0):return np.nan
    return np.trapezoid(q[col],q.time)/(b-a)
def minimum(d,c):
    return (np.nan,np.nan) if d.empty else (d[c].min(),d.time[d[c].idxmin()])
def csv(name,records):
    pd.DataFrame(records).to_csv(OUT/(name+'.csv'),index=False,encoding='utf-8-sig')
def independent(d,m,case,arm):
    t=m['Ts'];x=m['Xs'];R=m['C']['R'];A=m['C']['A'];P=R['P'];O=m['O']
    assert len(d)==len(x) and x.shape[1]==36 and np.max(abs(t-d.time))<1e-10
    assert np.all(np.diff(t)>0) and np.isfinite(x).all()
    rot=Rotation.from_quat(x[:,[7,8,9,6]])
    vel=rot.apply(x[:,3:6]);vg=np.hypot(vel[:,0],vel[:,1]);v=np.linalg.norm(vel,axis=1)
    gam=np.arctan2(-vel[:,2],vg);chi=np.arctan2(vel[:,1],vel[:,0])
    xf=x[:,19:25];eg=gam-xf[:,1];ec=(chi-xf[:,2]+np.pi)%(2*np.pi)-np.pi;ev=v-xf[:,0]
    q=np.clip((t-3)/10,0,1);lo={'C1':31,'C2':31.5}[case];v0=R['Vb']
    vreq=v0+(lo-v0)*ramp(q);integ=10*(2.5*q**4-3*q**5+q**6)+np.maximum(0,t-13)
    rd=np.column_stack([v0*t+(lo-v0)*integ,np.zeros_like(t),np.full_like(t,-R['h_ref'])])
    pos=np.linalg.norm(x[:,:3]-rd,axis=1)
    vref=xf[:,0]*np.cos(xf[:,1])
    vrefg=np.hypot(vref*np.cos(xf[:,2])+x[:,16],vref*np.sin(xf[:,2])+x[:,17])
    assert np.max(abs(x[:,16:19]))<1e-10
    z=np.clip((t-3)/22,0,1);vcenter=v0+(31.25-v0)*ramp(z)
    c=np.column_stack([vcenter,np.zeros_like(t),(31.25-v0)*30*z**2*(1-z)**2/22,np.zeros((len(t),2))])
    u=xf[:,[0,1,3,4,5]]-c;h=.01-np.einsum('ij,jk,ik->i',u,R['PB'],u)
    if arm=='Direct':assert np.max(abs(u))<1e-12 and np.all(xf[:,2]==0)
    def force(vv,gg,vr,gr,cr):
        mass=P['m']; X=mass*(vr+R['g']*np.sin(gg));N=mass*np.hypot(vv*gr+R['g']*np.cos(gg),vv*np.cos(gg)*cr)
        def mapped(a):
            k=.5*R['rho']*P['S']*vv**2;th=P['theta']
            cd=np.polynomial.polynomial.polyval(a,th[[0,1,4,7,9]])
            cl=np.polynomial.polynomial.polyval(a,th[[15,16,20,21,22]])
            return k*cl+(k*cd+X)*np.tan(a)
        return N-mapped(R['alphaLo']),mapped(R['alphaHi'])-N
    flo,fhi=force(v,gam,xf[:,3]-A['kV']*ev,xf[:,4]-A['kGamma']*eg,xf[:,5]-A['kChi']*ec)
    rlo,rhi=force(xf[:,0],xf[:,1],*xf[:,3:6].T)
    reconstructed=dict(vground=vg,vrequest=vreq,vref_ground=vrefg,positionError=pos,airSpeed=v,
        eV=ev,egamma_deg=np.rad2deg(eg),echi_deg=np.rad2deg(ec),refLo=rlo,refHi=rhi,fbLo=flo,fbHi=fhi,Tact=x[:,35])
    residuals={k:float(np.max(abs(d[k]-out))) for k,out in reconstructed.items()}
    assert max(residuals.values())<1e-8,residuals
    if arm=='F':assert np.max(abs(d.hmove-h))<1e-10
    d['tube_h']=h;d['err']=abs(vg-vreq);d['ref_err']=abs(vrefg-vreq)
    return residuals
def describe(d,case,arm,step,status,coverage):
    full=abs(d.time.iloc[-1]-35)<1e-9 and status=='completed'
    crit=d[d.time.between(14,20)]
    r=dict(case=case,arm=arm,step=step,status=status,end=d.time.iloc[-1],samples=len(d),
        actual_MAE=averaged(d,'err',3,16),reference_MAE=averaged(d,'ref_err',3,16),
        position_RMS=np.sqrt(np.mean(d.positionError**2)) if full else np.nan,position_peak=d.positionError.max(),
        plateau_actual=averaged(d,'vground',25,35),plateau_reference=averaged(d,'vref_ground',25,35),
        chartVmin=d.chartV.min(),chartVhmin=d.chartVh.min())
    r['plateau_bias']=r['plateau_actual']-{'C1':31,'C2':31.5}[case]
    for suffix,window in [('critical',crit),('observed',d)]:
        for col in ['fbHi','fbLo']:
            r[f'{col}_{suffix}'],r[f'{col}_{suffix}_t']=minimum(window,col)
            if suffix=='critical' and (window.empty or window.time.iloc[-1]<20):
                r[f'{col}_{suffix}']=r[f'{col}_{suffix}_t']=np.nan
        if window.empty:continue
        tube=window.tube_h>=0
        box=(abs(window.eV)<=.125)&(abs(window.egamma_deg)<=.45)&(abs(window.echi_deg)<=.35)
        eg=tube&box;joint=eg&(window.theta_sum_deg<=.70)
        coverage.append(dict(case=case,arm=arm,step=step,window=suffix,n=len(window),tube=tube.sum(),box=box.sum(),EG=eg.sum(),J=joint.sum(),min_tube_h=window.tube_h.min(),
            max_eV=abs(window.eV).max(),max_eg_deg=abs(window.egamma_deg).max(),max_ec_deg=abs(window.echi_deg).max(),max_theta_sum_deg=window.theta_sum_deg.max()))
        if suffix=='critical':
            assert (window.loc[eg,['fbHi','fbLo','normalMargin','bankMargin']]>0).all().all() and window.Traw[joint].between(0,60,inclusive='neither').all()
    for col in ['Tref','Traw','Tact','alphaActual_deg']:
        r[col+'_min']=d[col].min();r[col+'_max']=d[col].max()
    for col in ['surfaceAmp','surfaceRate','thrustAmp','thrustRate','surfaceCmdAtBound','thrustCmdClip']:
        r[col+'_samples']=int(d[col].sum())
    return r

def comparison(label,key,x,y,sign,units,floor=None,scale=None):
    z=sign*(x[key]-y[key]);gate=None;resolved=None;large=None
    if floor is not None:
        gate=5*max(abs(x.loc['005',key]-x.loc['010',key]),abs(y.loc['005',key]-y.loc['010',key]),floor)
        if not np.isfinite(gate):gate=None
        if gate is not None and np.isfinite(z).all():
            resolved=bool(z.prod()>0 and z.abs().min()>gate)
            large=bool(resolved and z.abs().min()>=scale) if scale is not None else None
    return dict(contrast=label,metric=key,units=units,coarse=z['010'],fine=z['005'],drift_gate=gate,resolved=resolved,absolute_research_scale=scale,exceeds_scale=large)
def run():
    protocol=json.loads((H/'protocol.json').read_bytes())
    for p,digest in protocol['hashes'].items():assert sha(H.parent/p)==digest,p
    metrics=[];cover=[];resid=[];traces={};states={};events=[];common=[]
    for sc in ['C1','C2']:
        for st in ['010','005']:
            for arm in ['Direct','F']:
                folder=H/'runs'/f'{sc}_{st}' if arm=='Direct' else PREV/'runs'/f'{sc}_F_{st}'
                m=loadmat(folder/('direct.mat' if arm=='Direct' else 'F.mat'),simplify_cells=True)
                d=pd.read_csv(folder/('trace.csv' if arm=='Direct' else 'replay.csv'),float_precision='round_trip')
                assert m['C']['MaxStep']==int(st)/1000
                assert m['O']['SpeedStart']==3 and m['O']['SpeedRise']==10 and m['O']['Vlow']==protocol['cases'][sc]
                rr=independent(d,m,sc,arm);resid.append(dict(case=sc,arm=arm,step=st,**rr))
                status=m['Summary']['status'];metrics.append(describe(d,sc,arm,st,status,cover))
                traces[sc,st,arm]=d;states[sc,st,arm]=m['Xs']
                indices=m['Summary']['event_indices']
                times=m['Summary']['events_s' if arm=='Direct' else 'event_times_s']
                for i,t in zip(np.atleast_1d(indices),np.atleast_1d(times)):
                    events.append(dict(case=sc,step=st,arm=arm,index=int(i),time=float(t)))
            a,b=[traces[sc,st,arm] for arm in ['Direct','F']]
            x=a.merge(b,on='time',suffixes=('_direct','_F'),validate='one_to_one')
            common.append(dict(case=sc,step=st,from_s=x.time.min(),to_s=x.time.max(),n=len(x),
                direct_position_RMS=np.sqrt(np.mean(x.positionError_direct**2)),F_position_RMS=np.sqrt(np.mean(x.positionError_F**2))))
    tab=pd.DataFrame(metrics).set_index(['case','arm','step']).sort_index()
    prev=pd.read_csv(PREV/'results/run_metrics.csv',dtype={'step':str}).set_index(['case','arm','step']).sort_index()
    contrasts=[];shares=[];qa=[]
    for sc in ['C1','C2']:
        direct=tab.loc[sc,'Direct'];f=tab.loc[sc,'F']
        for col,old in [('actual_MAE','MAE_3_16'),('plateau_actual','vground_25_35'),('position_RMS','position_RMS35'),('fbHi_critical','min_fb_hi_window'),('fbLo_critical','min_fb_lo_window')]:
            assert (abs(f[col]-prev.loc[sc,'F'][old])<1e-8).all()
        for key,units,floor,scale in [('actual_MAE','m/s',1e-4,.606),('reference_MAE','m/s',1e-4,.606),('position_RMS','m',None,None),('position_peak','m',None,None),
            ('fbHi_critical','N',1e-7,1),('fbLo_critical','N',1e-7,1),('fbHi_observed','N',1e-7,1),('fbLo_observed','N',1e-7,1),('plateau_actual','m/s',1e-4,None)]:
            contrasts.append(dict(case=sc,**comparison('Direct_minus_F',key,direct,f,1,units,floor,scale)))
        for st in ['010','005']:
            od=prev.loc[(sc,'O',st),'MAE_3_16'];v=direct.loc[st,'actual_MAE'];w=f.loc[st,'actual_MAE']
            shares.append(dict(case=sc,step=st,O=od,Direct=v,F=w,O_minus_F=od-w,O_minus_Direct=od-v,Direct_minus_F=v-w,
                              fraction_of_O_F_gain_at_Direct=(od-v)/(od-w),fraction_remaining=(v-w)/(od-w)))
    for arm in ['Direct','F']:
        x=tab.xs(('C2',arm));y=tab.xs(('C1',arm))
        for key in ['plateau_actual','plateau_reference']:
            contrasts.append(dict(case=arm,**comparison('C2_minus_C1',key,x,y,1,'m/s',1e-4,None)))
    for st in ['010','005']:
        a,b=[states[sc,st,'Direct'] for sc in ['C1','C2']]
        t1,t2=[traces[sc,st,'Direct'].time for sc in ['C1','C2']]
        _,i,j=np.intersect1d(t1,t2,return_indices=True)
        # Exclude the passive position guidance integral; the physical plant receives the same analytic reference.
        body=np.max(abs(a[i,:13]-b[j,:13]));servo=np.max(abs(a[i,32:36]-b[j,32:36]))
        qa.append(dict(step=st,n=len(i),max_abs_13plant_state_diff=body,max_abs_4actuator_state_diff=servo,
            max_abs_passive_integral_diff=np.max(abs(a[i,13:16]-b[j,13:16]))))
    for name,records in [('run_metrics',metrics),('paired_differences',contrasts),('speed_comparison',shares),('coverage',cover),('common_windows',common),('reconstruction_residuals',resid),('cross_request_checks',qa)]:
        csv(name,records)
    pd.DataFrame(events,columns=['case','step','arm','index','time']).to_csv(OUT/'events.csv',index=False)
    return tab,traces
if __name__=='__main__':run()
