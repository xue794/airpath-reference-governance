"""Assemble the nine fixed cells without integrating or changing the WP3 evaluator."""
from pathlib import Path
import json,csv,sys,hashlib
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H/'analysis'))
from interval_wp3 import Evaluator

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item()),encoding='utf-8')
def csvread(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def writecsv(p,x):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(x[0]));w.writeheader();w.writerows(x)
def segments(t,m):
    out=[];start=None
    for k,yes in enumerate(m):
        if yes and start is None:start=k
        if start is not None and (not yes or k==len(m)-1):
            last=k if yes else k-1
            out.append(dict(first_sample_s=float(t[start]),last_sample_s=float(t[last]),samples=last-start+1));start=None
    return out

new=read(H/'analysis/results/transfer_summary.json');old=read(H/'analysis/prior_results/transfer_summary.json')
transfer={x['key']:x for x in old+new}
layered={x['Method']+'_'+x['Case']:x for x in read(H/'results/case_summary.json')}
assert len(transfer)==len(layered)==9
table=[];series={};uncovered={};fa=read(H/'analysis/prior_results/joint_failure_attribution.json')
oldsummary={x['method']+'_'+x['case']:x for x in read(H.parent/'revision_wp3_20261107/results/summary_table.json')}
diffs=[]
for case in ['hold','varying','noisy']:
    for method in ['Transparent','Proposed','MatchedQP']:
        key=method+'_'+case;s=transfer[key];z=layered[key]
        folder=H/('analysis/results' if method=='Transparent' and case!='varying' else 'analysis/prior_results')
        rows=csvread(folder/f'{key}_transfer.csv')
        L=loadmat(H/'results'/f'{key}_layered.mat',simplify_cells=True)
        v={str(n):L['Values'][:,i] for i,n in enumerate(L['valueNames'])}
        d={n:np.array([float(x[n]) for x in rows]) for n in ['time_s','error_box_mps','radius_J_mps','theta_sum_deg']}
        covered=np.array([x['inclusion_J_S4']=='True' for x in rows]);t=d['time_s']
        assert len(t)==z['Samples']==s['samples'] and np.max(abs(t-L['t']))<1e-12
        assert np.max(abs(np.array([float(x['Traw_N']) for x in rows])-v['raw_T_N']))<1e-10
        series[key]=(v,d,covered);uncovered[key]=segments(t,~covered)
        rr=s['radius_stats']['J']
        row=dict(case=case,method=method,new_integration=method=='Transparent' and case!='varying',samples=s['samples'],completed=z['Completed'],end_s=s['end_s'],
            RMS_grid_m=z['PositionRMS_grid_m'],RMS_time_m=z['PositionRMS_time_m'],RMS_late_m=z['PositionRMS_late_m'],peak_error_m=z['PositionPeak_m'],
            reference_T_min_N=z['ReferenceThrustMin_N'],reference_T_max_N=z['ReferenceThrustMax_N'],
            reference_thrust_reserve_violations=z['ReferenceThrustReserveViolationSamples'],reference_endpoint_reserve_violations=z['ReferenceEndpointReserveViolationSamples'],
            reference_other_state_violations=int(np.any(L['Ref'][:,0:23]<-1e-8,axis=1).sum()),
            derivative_nonbackup_violations=z['RawDerivativeNonbackupViolationSamples'],normalized_derivative_violations=z['NormalizedDerivativeViolationSamples'],
            backup_outside_samples=z['BackupOutsideSamples'],backup_derivative_violations=z['BackupDerivativeViolationSamples'],max_VB_over_rho=z['MaxVBoverRho'],
            raw_T_min_N=z['RawThrustMin_N'],raw_T_max_N=z['RawThrustMax_N'],raw_T_violations=z['FeedbackRawThrustViolationSamples'],
            feedback_margin_violations=int(np.any(L['Feed'][:,[0,1,2,3,4,11,12]]<-1e-8,axis=1).sum()),
            feedback_reserve_violations=int(np.any(L['Feed'][:,[5,6,13,14]]<-1e-8,axis=1).sum()),
            modeled_actuator_violations=z['ModeledActuatorViolationSamples'],study_domain_violations=z['StudyDomainViolationSamples'],
            thrust_target_clipped=z['ThrustTargetClippedSamples'],thrust_rate_limited=z['ThrustRateLimitedSamples'],surface_rate_limited=z['SurfaceRateLimitedSamples'],allocator_active=z['AllocatorBoundActiveSamples'],
            max_allocation_residual=float(v['allocation_residual_norm'].max()),max_servo_residual=float(v['servo_residual_norm'].max()),
            J_min_mps=rr['min'],J_median_mps=rr['median'],J_max_mps=rr['max'],J_positive=s['counts']['certified_J'],
            J_inclusion=int(covered.sum()),J_percent=100*float(covered.mean()),G_inclusion=s['counts']['inclusion_G'],E_inclusion=s['counts']['inclusion_E'],
            max_error_box_mps=s['max_error_box_mps'],max_theta_sum_deg=s['max_theta_sum_deg'],hold_max_state_difference=z['HoldStateDifference'],
            min_reference_alpha_deg=float(v['reference_alpha_deg'].min()),max_reference_alpha_deg=float(v['reference_alpha_deg'].max()),
            min_feedback_alpha_deg=float(v['command_alpha_deg'].min()),max_feedback_alpha_deg=float(v['command_alpha_deg'].max()),
            min_actual_alpha_deg=float(v['actual_alpha_deg'].min()),max_actual_alpha_deg=float(v['actual_alpha_deg'].max()))
        table.append(row)
        if key in oldsummary:
            prev=oldsummary[key]
            delta=max(abs(row['RMS_grid_m']-prev['RMS_grid_m']),abs(row['raw_T_min_N']-prev['observed_min_Traw_N']),abs(row['reference_T_min_N']-prev['observed_min_reference_T_N']))
            assert delta<1e-10 and row['J_inclusion']==prev['joint_S4_count']
            diffs.append(dict(key=key,max_metric_difference=delta))
        else:
            X=loadmat(H/'analysis/results'/f'{key}_export.mat',simplify_cells=True);E=Evaluator(X['C']);D=X['D'];ii={str(n):i for i,n in enumerate(X['names'])};fails=[]
            for k in np.flatnonzero(~covered):
                E.prepare(D[k,1:7],D[k,ii['ref_alpha']]);b=E.bounds(d['error_box_mps'][k],thrust=True)
                fails.append(dict(time_s=t[k],raw_thrust_lower=b['thrust_lower'],raw_thrust_upper=b['thrust_upper'],endpoint_pass=bool(b['E']),geometry_pass=bool(b['G']),
                    lower_fails=b['thrust_lower'] is None or b['thrust_lower']<=0,upper_fails=b['thrust_upper'] is None or b['thrust_upper']>=60,kappa_fails=b['kappa_lower'] is None or b['kappa_lower']<=0))
            if fails:writecsv(H/'results'/f'{key}_uncovered_bounds.csv',fails)
            fa.append(dict(key=key,uncovered_samples=int((~covered).sum()),lower_thrust_failures=sum(x['lower_fails'] for x in fails),upper_thrust_failures=sum(x['upper_fails'] for x in fails),
                geometry_failures=sum(not x['geometry_pass'] for x in fails),kappa_failures=sum(x['kappa_fails'] for x in fails)))
dump(H/'results/summary_table.json',table);writecsv(H/'results/summary_table.csv',table)
dump(H/'results/combined_transfer_summary.json',list(transfer.values()));dump(H/'results/joint_failure_attribution.json',fa);dump(H/'results/uncovered_sample_segments.json',uncovered)
dump(H/'checks/assembly_check.json',dict(cells=len(table),new_integrations=2,reused_cells=7,samples=sum(x['samples'] for x in table),
    old_seven_metric_checks=diffs,all_complete=all(x['completed'] for x in table),all_EG_included=all(x['E_inclusion']==x['G_inclusion']==x['samples'] for x in table),
    all_J_positive=all(x['J_positive']==x['samples'] for x in table),all_theta_condition_observed=all(x['max_theta_sum_deg']<=.7 for x in table),
    max_reconstruction_difference=max(max(x['reconstruction_max_errors']) for x in transfer.values())))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
colors={'Transparent':'#777777','Proposed':'#2166ac','MatchedQP':'#b35806'}
labels={'Transparent':'Transparent','Proposed':'Scalar','MatchedQP':'Matched-QP'}
fig,axs=plt.subplots(3,1,figsize=(9,7),sharex=True,layout='constrained')
for method in colors:
    v,d,c=series[method+'_noisy'];t=d['time_s']
    for a,col in zip(axs,['reference_T_N','raw_T_N','actuator_T_N']):a.plot(t,v[col],color=colors[method],lw=1.1,label=labels[method])
for a,title in zip(axs,['Reference thrust','Raw feedback thrust before actuator limiting','Modeled actuator thrust']):
    for z in [0,60]:a.axhline(z,color='#222222',ls='--',lw=.7)
    a.set_ylim(-2,62);a.set_ylabel('Thrust (N)');a.set_title(title,fontsize=10);a.grid(alpha=.18)
for z in [3,57]:axs[0].axhline(z,color='#a83b38',ls=':',lw=1)
axs[0].legend(ncol=3,fontsize=9,loc='upper right');axs[-1].set_xlabel('Time (s)')
fig.savefig(H/'figures/noisy_thrust_layers.png',dpi=200);fig.savefig(H/'figures/noisy_thrust_layers.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(10,3.4),layout='constrained')
for method in colors:
    v,d,c=series[method+'_noisy'];t=d['time_s']
    axs[0].plot(t,v['position_error_m'],color=colors[method],lw=1,label=labels[method])
    axs[1].plot(t,v['VB_over_rho'],color=colors[method],lw=1)
axs[0].set_ylabel('Position-error norm (m)');axs[0].legend(fontsize=8)
axs[1].set_ylabel('Backup value / rho');axs[1].set_yscale('log');axs[1].axhline(1,color='black',ls='--',lw=.8)
for a in axs:a.set_xlabel('Time (s)');a.grid(alpha=.15)
fig.savefig(H/'figures/noisy_tracking_backup.png',dpi=210);fig.savefig(H/'figures/noisy_tracking_backup.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(10,4.1),sharey=True,layout='constrained')
for a,case in zip(axs,['hold','noisy']):
    v,d,c=series['Transparent_'+case];t=d['time_s']
    radii=d['radius_J_mps'].copy();uncert=radii<=0;radii[uncert]=np.nan
    a.plot(t,radii,color='#2166ac',label='Positive joint sufficient half-width')
    a.plot(t,d['error_box_mps'],color='#222222',lw=1,label='Observed scaled error')
    a.fill_between(t,0,1,where=~c,color='#c85d44',alpha=.13,step='mid',label='Joint condition not established')
    if np.any(uncert):a.scatter(t[uncert],np.zeros(int(uncert.sum())),color='#a83b38',s=13,marker='x',clip_on=False,label='No positive certificate')
    a.set_title('Transparent '+case);a.set_ylim(0,1);a.set_xlabel('Time (s)');a.grid(alpha=.15)
axs[0].set_ylabel('Half-width / scaled error (m/s)')
handles,labs=axs[1].get_legend_handles_labels();fig.legend(handles,labs,fontsize=8,loc='outside upper center',ncol=2)
fig.savefig(H/'figures/transparent_transfer.png',dpi=210);fig.savefig(H/'figures/transparent_transfer.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(10,3.4),sharex=True,layout='constrained')
v,d,c=series['Transparent_noisy'];t=d['time_s'];mask=(t>=23)&(t<=31)
axs[0].plot(t[mask],v['reference_T_N'][mask],color='#777777',lw=1.4,label='Reference thrust')
axs[0].axhline(3,color='#a83b38',ls=':',label='3 N reference reserve');axs[0].axhline(0,color='black',ls='--',label='0 N reference bound')
axs[0].scatter(t[v['reference_T_N']<0],v['reference_T_N'][v['reference_T_N']<0],s=10,color='#b23e32',label='17 negative saved points')
axs[1].plot(t[mask],v['raw_T_N'][mask],color='#2166ac',label='Raw feedback thrust')
axs[1].plot(t[mask],v['actuator_T_N'][mask],color='#b35806',ls='--',label='Actuator thrust')
axs[0].set_ylim(-1,7);axs[1].set_ylim(0,25)
for a in axs:a.set_xlabel('Time (s)');a.set_ylabel('Thrust (N)');a.set_xlim(23,31);a.legend(fontsize=8);a.grid(alpha=.15)
fig.savefig(H/'figures/noisy_negative_reference.png',dpi=210);fig.savefig(H/'figures/noisy_negative_reference.svg');plt.close(fig)
print(json.dumps(dict(cells=len(table),samples=sum(x['samples'] for x in table),new=[x for x in table if x['new_integration']]),ensure_ascii=False))
