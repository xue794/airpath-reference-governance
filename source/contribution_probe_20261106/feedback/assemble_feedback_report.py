"""Assemble tables and a standalone scientific figure from the audit outputs."""
from pathlib import Path
import csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]/'endpoint_probe_20261106'/'full6dof'
def table(path):
    with path.open(newline='',encoding='utf-8') as f:
        return [{k:(v=='True' if v in ['True','False'] else float(v)) for k,v in r.items()} for r in csv.DictReader(f)]
summaries=json.loads((ROOT/'transfer_summary.json').read_text())
post=json.loads((ROOT/'posthoc_attitude_summary.json').read_text())
combined=[]
for s in summaries:
    c=s['case'];m=s['method'];r=table(ROOT/f'{c}_{m}_transfer.csv');n=len(r)
    V=lambda k:np.array([x[k] for x in r])
    L=loadmat(SOURCE/c/(m+'.mat'),simplify_cells=True)
    row=dict(case=c,method=m,samples=n,end_time_s=s['end_time_s'],
        min_reference_lower_N=V('reference_lower_N').min(),min_reference_upper_N=V('reference_upper_N').min(),
        min_feedback_lower_N=V('feedback_lower_N').min(),min_feedback_upper_N=V('feedback_upper_N').min(),
        min_reference_thrust_N=L['D'][:,9].min(),max_reference_thrust_N=L['D'][:,9].max(),
        min_Tstar_N=V('Tstar_N').min(),max_Tstar_N=V('Tstar_N').max(),
        min_Traw_N=V('Traw_N').min(),max_Traw_N=V('Traw_N').max(),
        min_endpoint_radius_mps=V('radius_endpoints_mps').min(),
        min_endpoint_joint_angle_bound_deg=V('radius_endpoints_mps').min()/38.1*180/math.pi,
        endpoint_radius_capped_at_1_samples=int(sum(V('radius_endpoints_mps')==1)),
        endpoint_inclusion_samples=int(sum(V('endpoint_tube_inclusion'))),endpoint_inclusion_pct=100*sum(V('endpoint_tube_inclusion'))/n,
        min_physical_theta0p7_radius_mps=V('radius_four_physical_mps').min(),
        physical_theta0p7_inclusion_samples=int(sum(V('physical_tube_inclusion'))),
        physical_theta0p7_inclusion_pct=100*sum(V('physical_tube_inclusion'))/n,
        min_reserved_theta0p7_radius_mps=V('radius_four_reserved_mps').min(),
        reserved_theta0p7_inclusion_samples=int(sum(V('reserved_tube_inclusion'))),
        attitude_error_le0p7_samples=int(sum(V('attitude_error_deg')<=.7)),
        max_abs_eV_mps=max(abs(V('eV_mps'))),max_abs_egamma_deg=max(abs(V('egamma_deg'))),max_abs_echi_deg=max(abs(V('echi_deg'))),
        max_attitude_error_deg=max(V('attitude_error_deg')),
        max_error_box_mps=max(V('error_box_mps')),max_upper_margin_loss_N=max(V('loss_upper_N')),
        min_four_physical_theta2p5_radius_mps=None,physical_theta2p5_inclusion_samples=None,physical_theta2p5_inclusion_pct=None)
    if m=='N10':
        p=next(p for p in post if p['case']==c)
        row.update(min_four_physical_theta2p5_radius_mps=p['minimum_radius_mps'],physical_theta2p5_inclusion_samples=p['inclusion_samples'],physical_theta2p5_inclusion_pct=100*p['inclusion_samples']/n)
    combined.append(row)
with (ROOT/'comparison_summary.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=combined[0]);w.writeheader();w.writerows(combined)
(ROOT/'comparison_summary.json').write_text(json.dumps(combined,indent=2),encoding='utf-8')
fig,ax=plt.subplots(3,2,figsize=(12,9),sharex='col',layout='constrained')
for j,c in enumerate(['V31_C60_warm_fine','V32_C60_warm_fine']):
    r=table(ROOT/f'{c}_N10_transfer.csv');p=table(ROOT/f'{c}_N10_theta2p5.csv')
    V=lambda k:np.array([x[k] for x in r]);t=V('time_s');D=loadmat(SOURCE/c/'N10.mat',simplify_cells=True)['D']
    ax[0,j].plot(t,V('reference_upper_N'),label='Reference upper endpoint',color='#4676b7')
    ax[0,j].plot(t,V('feedback_upper_N'),label='Feedback upper endpoint',color='#bd4528')
    ax[0,j].axhline(0,color='.4',lw=.8);ax[0,j].axhline(10,color='.6',ls=':',lw=.8)
    ax[0,j].set_ylim(-1,30);ax[0,j].set_title(f'{31+j} m/s request; N10');ax[0,j].set_ylabel('Upper endpoint margin (N)')
    ax[1,j].plot(t,V('error_box_mps'),label='Observed error scale',color='black',lw=1.5)
    ax[1,j].plot(t,V('radius_endpoints_mps'),label='Endpoint sufficient radius',color='#4676b7')
    ax[1,j].plot(t,[x['radius_four_physical_theta2p5_mps'] for x in p],label='Physical four-margin radius\n(post hoc attitude <= 2.5 deg)',color='#548843')
    ax[1,j].set_ylim(0,.62);ax[1,j].set_ylabel('Weighted error scale (m/s)')
    ax[2,j].plot(t,D[:,9],label='Reference thrust',color='#4676b7')
    ax[2,j].plot(t,V('Traw_N'),label='Feedback raw thrust',color='#bd4528')
    ax[2,j].axhline(0,color='.3',lw=.8);ax[2,j].axhline(60,color='.3',lw=.8)
    ax[2,j].axhline(3,color='.6',ls=':',lw=.8);ax[2,j].axhline(57,color='.6',ls=':',lw=.8)
    ax[2,j].set_ylim(-1,62);ax[2,j].set_ylabel('Thrust (N)');ax[2,j].set_xlabel('Time (s)')
    for k in range(3):ax[k,j].grid(alpha=.18);ax[k,j].set_xlim(0,35)
for k in range(3):ax[k,0].legend(fontsize=8,loc='upper right')
fig.suptitle('Conditional reference-to-feedback transfer audit (archived fine trajectories)',fontsize=13)
fig.savefig(ROOT/'feedback_transfer_audit.png',dpi=180)
fig.savefig(ROOT/'feedback_transfer_audit.pdf')
print(json.dumps(combined,indent=2))
