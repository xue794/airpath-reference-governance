"""Summaries, failure attribution and scientific figures from the fixed audit."""
from pathlib import Path
import json,csv,math,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from interval_wp3 import Evaluator
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item()),encoding='utf-8')
def writecsv(p,rows):
    with p.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
def csvread(p):
    with p.open(encoding='utf-8') as f:return list(csv.DictReader(f))
def segments(t,mask):
    groups=[];first=None
    for i,yes in enumerate(mask):
        if yes and first is None:first=i
        if first is not None and (not yes or i==len(mask)-1):
            last=i if yes else i-1;groups.append(dict(first_sample_s=float(t[first]),last_sample_s=float(t[last]),sample_count=last-first+1));first=None
    return groups
summaries=json.loads((ROOT/'results/transfer_summary.json').read_text(encoding='utf-8'))
export=json.loads((ROOT/'results/export_summary.json').read_text(encoding='utf-8'))
ex={v['Method']+'_'+v['Case']:v for v in export}
table=[];margins=[];attribution=[];series={};allseg={};minimum_details=[]
for s in summaries:
    key=s['key'];rows=csvread(ROOT/'results'/f'{key}_transfer.csv');L=loadmat(ROOT/'results'/f'{key}_export.mat',simplify_cells=True)
    D=L['D'];names=[str(x) for x in L['names']];d={n:D[:,i] for i,n in enumerate(names)};E=Evaluator(L['C']);z=ex[key]
    vals={k:np.array([float(x[k]) for x in rows]) for k in ['time_s','radius_E_mps','radius_G_mps','radius_J_mps','error_box_mps','theta_sum_deg','theta_direct_deg']}
    covered=np.array([x['inclusion_J_S4']=='True' for x in rows]);t=vals['time_s'];series[key]=(vals,d,covered)
    allseg[key]=segments(t,~covered)
    failures=[]
    for k in np.flatnonzero(~covered):
        E.prepare(D[k,1:7],d['ref_alpha'][k]);b=E.bounds(vals['error_box_mps'][k],thrust=True)
        failures.append(dict(time_s=t[k],raw_thrust_lower=b['thrust_lower'],raw_thrust_upper=b['thrust_upper'],
            endpoint_pass=bool(b['E']),geometry_pass=bool(b['G']),lower_fails=b['thrust_lower'] is None or b['thrust_lower']<=0,
            upper_fails=b['thrust_upper'] is None or b['thrust_upper']>=60,kappa_fails=b['kappa_lower'] is None or b['kappa_lower']<=0))
    if failures:writecsv(ROOT/'results'/f'{key}_uncovered_bounds.csv',failures)
    fa=dict(key=key,uncovered_samples=int((~covered).sum()),lower_thrust_failures=sum(x['lower_fails'] for x in failures),
        upper_thrust_failures=sum(x['upper_fails'] for x in failures),geometry_failures=sum(not x['geometry_pass'] for x in failures),
        kappa_failures=sum(x['kappa_fails'] for x in failures))
    attribution.append(fa)
    n=s['samples'];r=s['radius_stats']['J'];late=t>=(30 if s['case']=='hold' else 60)
    table.append(dict(case=s['case'],method=s['method'],samples=n,J_min_mps=r['min'],J_median_mps=r['median'],J_max_mps=r['max'],
        joint_S4_count=int(covered.sum()),joint_S4_percent=100*covered.mean(),G_count=s['counts']['inclusion_G'],
        late_joint_count=int(covered[late].sum()),late_samples=int(late.sum()),
        max_error_box_mps=s['max_error_box_mps'],max_theta_sum_deg=s['max_theta_sum_deg'],
        RMS_grid_m=z['PositionRMS'],RMS_time_m=float(np.sqrt(np.trapezoid(d['position_error']**2,t)/(t[-1]-t[0]))),
        observed_min_Traw_N=float(d['Traw'].min()),observed_max_Traw_N=float(d['Traw'].max()),
        observed_min_reference_T_N=float(d['ref_T'].min()),observed_max_reference_T_N=float(d['ref_T'].max()),
        reference_thrust_reserve_violations=int(((d['ref_T']<3-1e-8)|(d['ref_T']>57+1e-8)).sum()),
        feedback_thrust_reserve_violations=int(((d['Traw']<3-1e-8)|(d['Traw']>57+1e-8)).sum()),
        backup_outside_samples=int((d['ref_hB']<-1e-8).sum()),
        modeled_actuator_violation_samples=z['PhysicalActuatorViolationSamples'],
        modeled_study_domain_violation_samples=z['StudyDomainViolationSamples'],
        input_source=z['Source']))
    mi=int(vals['radius_J_mps'].argmin());mx=int((vals['error_box_mps']/vals['radius_J_mps']).argmax())
    minimum_details.append(dict(key=key,minimum_radius_sample=dict(index=mi,time_s=t[mi],radius_mps=r['min'],reference_T=d['ref_T'][mi],raw_T=d['Traw'][mi],observed_error_mps=vals['error_box_mps'][mi]),
        maximum_relative_error_sample=dict(index=mx,time_s=t[mx],radius_mps=vals['radius_J_mps'][mx],observed_error_mps=vals['error_box_mps'][mx],ratio=vals['error_box_mps'][mx]/vals['radius_J_mps'][mx])))
    for category,M,nn,units in [('reference',L['Ref'],['normal','bank_plus','bank_minus','endpoint_lower','endpoint_upper','thrust_lower','thrust_upper'],['m/s2']*3+['N']*4),
        ('feedback',L['Feed'],['normal','bank_plus','bank_minus','endpoint_lower','endpoint_upper','thrust_lower','thrust_upper'],['m/s2']*3+['N']*4)]:
        for j,name in enumerate(nn):
            k=int(M[:,j].argmin());bad=M[:,j]<-1e-8;ix=np.flatnonzero(bad)
            margins.append(dict(case=s['case'],method=s['method'],layer=category,quantity=name,unit=units[j],minimum=float(M[k,j]),time_s=float(t[k]),negative_samples=int(bad.sum()),first_negative_s=None if not len(ix) else float(t[ix[0]]),tolerance=1e-8))
    print('ASSEMBLED',key,fa,flush=True)
dump(ROOT/'results/summary_table.json',table);writecsv(ROOT/'results/summary_table.csv',table)
dump(ROOT/'results/margin_summary.json',margins);writecsv(ROOT/'results/margin_summary.csv',margins)
dump(ROOT/'results/joint_failure_attribution.json',attribution);dump(ROOT/'results/uncovered_sample_segments.json',allseg)
dump(ROOT/'results/minimum_radius_details.json',minimum_details)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','font.family':'DejaVu Sans'})
fig,axs=plt.subplots(2,3,figsize=(11.2,6.1),sharex=True,sharey='row',layout='constrained')
for col,(method,label) in enumerate([('Transparent','Transparent'),('Proposed','Scalar'),('MatchedQP','Matched-QP')]):
    v,d,covered=series[method+'_varying'];t=v['time_s'];a=axs[0,col]
    a.plot(t,v['radius_G_mps'],color='#89939b',ls='--',label='Geometry box (search cap)')
    a.plot(t,v['radius_J_mps'],color='#2166ac',label='Joint sufficient half-width')
    a.plot(t,v['error_box_mps'],color='#222222',lw=1,label='Observed scaled air-path error')
    a.fill_between(t,0,1.05,where=~covered,color='#c85d44',alpha=.1,step='mid',label='Joint condition not established')
    a.set_ylim(0,1.05);a.set_title(label);a.grid(alpha=.15)
    a=axs[1,col];a.plot(t,v['theta_sum_deg'],color='#222222',lw=1,label='Observed S4 attitude sum')
    a.plot(t,v['theta_direct_deg'],color='#2166ac',ls=':',lw=1,label='Observed direct attitude error')
    a.axhline(.7,color='#ab3b36',ls='--',label='Conditional bound 0.70 deg');a.set_ylim(0,.76);a.grid(alpha=.15);a.set_xlabel('Time (s)')
axs[0,0].set_ylabel('Scaled half-width / error (m/s)');axs[1,0].set_ylabel('Attitude error (deg)')
handles,labels=axs[0,1].get_legend_handles_labels();fig.legend(handles,labels,fontsize=8,loc='outside upper center',ncol=2)
axs[1,1].legend(fontsize=7,loc='upper center')
fig.savefig(ROOT/'figures/varying_transfer.png',dpi=210);fig.savefig(ROOT/'figures/varying_transfer.svg');plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(6.7,4.8),sharex=True,layout='constrained')
v,d,covered=series['Proposed_varying'];t=v['time_s']
axs[0].plot(t,v['radius_G_mps'],color='#89939b',ls='--',label='Geometry box at search cap')
axs[0].plot(t,v['radius_J_mps'],color='#2166ac',label='Joint sufficient half-width')
axs[0].plot(t,v['error_box_mps'],color='#222222',lw=1,label='Observed air-path error')
axs[0].fill_between(t,0,1.05,where=~covered,color='#c85d44',alpha=.12,step='mid')
axs[0].set_ylim(0,1.05);axs[0].set_ylabel('Half-width / error (m/s)');axs[0].legend(fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2)
axs[1].plot(t,v['theta_sum_deg'],color='#222222',label='Observed S4 attitude sum')
axs[1].plot(t,v['theta_direct_deg'],color='#2166ac',ls=':',label='Observed direct attitude error')
axs[1].axhline(.7,color='#ab3b36',ls='--',label='Conditional bound 0.70 deg')
axs[1].set_ylim(0,.76);axs[1].set_ylabel('Attitude error (deg)');axs[1].set_xlabel('Time (s)');axs[1].legend(fontsize=8,loc='upper right')
for a in axs:a.grid(alpha=.15)
fig.savefig(ROOT/'figures/varying_scalar_transfer.png',dpi=220);fig.savefig(ROOT/'figures/varying_scalar_transfer.svg');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(10.3,3.5),layout='constrained')
for method,label,color in [('Transparent','Transparent','#777777'),('Proposed','Scalar','#2166ac'),('MatchedQP','Matched-QP','#b35806')]:
    v,d,c=series[method+'_varying'];t=v['time_s'];axs[0].plot(t,d['ref_T'],label=label,color=color,lw=1.1);axs[1].plot(t,d['Traw'],label=label,color=color,lw=1.1)
for a in axs:
    for z in [3,57]:a.axhline(z,color='#ab3b36',ls=':',lw=1)
    for z in [0,60]:a.axhline(z,color='black',ls='--',lw=.7)
    a.grid(alpha=.15);a.set_xlabel('Time (s)');a.set_ylabel('Thrust (N)');a.set_ylim(-2,62)
axs[0].set_title('Reference thrust');axs[1].set_title('Raw feedback thrust before actuator clipping');axs[0].legend(fontsize=8)
fig.savefig(ROOT/'figures/varying_thrust.png',dpi=210);fig.savefig(ROOT/'figures/varying_thrust.svg');plt.close(fig)
dump(ROOT/'checks/assembly_check.json',dict(samples=sum(s['samples'] for s in summaries),all_G_covered=all(s['counts']['inclusion_G']==s['samples'] for s in summaries),
    all_radii_positive=all(s['counts']['certified_J']==s['samples'] for s in summaries),all_attitude_condition_observed=all(s['counts']['attitude_S4_ok']==s['samples'] for s in summaries),
    no_feedback_margin_violation=all(not any(s['FeedbackViolationCounts']) for s in export),integrations=0))
print('WP3_ASSEMBLY_COMPLETE',flush=True)
