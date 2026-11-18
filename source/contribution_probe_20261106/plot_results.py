"""Scientific figures for matched AoA comparison and guarantee coverage."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H=Path(__file__).resolve().parent
OLD=H.parent/'endpoint_probe_20261106'/'full6dof'
audit=json.loads((H/'comparison_audit.json').read_text())
assert audit['complete'] and audit['all_available_run_checks_pass']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.grid':True,'grid.alpha':.2,'savefig.facecolor':'white'})
def load(speed,method):
    base=OLD/f'V{speed}_C60_warm_fine' if method.startswith('N') else H/'full6dof'/f'V{speed}_C60_fine'
    return np.genfromtxt(base/f'{method}_samples.csv',delimiter=',',names=True)
style={'N0':('#d55e00','-'),'A0':('#333333','--'),'N10':('#0072b2','-'),'Am':('#009e73','--')}
labels={'N0':'Endpoint, 0 N','A0':'Fixed AoA, 0 N','N10':'Endpoint, 10 N','Am':'Moving AoA, matched 10 N'}
fig,axes=plt.subplots(2,2,figsize=(12,7.8),constrained_layout=True)
for col,speed in enumerate((31,32)):
    methods=('N0','A0','N10','Am') if speed==31 else ('N10','Am')
    for method in methods:
        a=load(speed,method);color,ls=style[method]
        for row,field in enumerate(('feedback_raw_upper_N','position_error_m')):
            ax=axes[row,col];ax.plot(a['time_s'],a[field],color=color,ls=ls,lw=1.7,label=labels[method])
            if a['time_s'][-1]<35:ax.plot(a['time_s'][-1],a[field][-1],'x',color=color,ms=8,mew=2)
    axes[0,col].set(title=f'{speed} m/s request: feedback upper margin',ylabel='N',xlim=(10,22),ylim=(-1,35))
    axes[1,col].set(title=f'{speed} m/s request: position error',ylabel='m',xlim=(0,35))
    axes[0,col].axhline(0,color='gray',lw=1)
    axes[0,col].legend(fontsize=8,loc='upper right')
    for row in (0,1):axes[row,col].set_xlabel('Time (s)')
fig.suptitle('Force-endpoint versus coupled AoA protection\nMatched 10 N sets give very similar closed-loop outcomes; no backup-set row',fontsize=13)
fig.savefig(H/'aoa_comparison.png',dpi=180);fig.savefig(H/'aoa_comparison.svg');plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
for speed,color in [(31,'#0072b2'),(32,'#d55e00')]:
    a=load(speed,'N10')
    axes[0].plot(a['time_s'],a['VB_over_rho'],label=f'{speed} m/s',color=color)
    tr=np.genfromtxt(H/'feedback'/f'V{speed}_C60_warm_fine_N10_transfer.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    r=np.asarray(tr['radius_endpoints_mps'],float);e=np.asarray(tr['error_box_mps'],float)
    ratio=np.divide(e,r,out=np.full_like(e,np.nan),where=r>0)
    axes[1].plot(tr['time_s'],ratio,label=f'{speed} m/s',color=color)
for ax in axes:
    ax.axhline(1,color='#333333',ls='--',lw=1);ax.set_xlabel('Time (s)');ax.set_yscale('log');ax.legend()
axes[0].set(title='Original certified backup-set coverage',ylabel='V_B / rho_B (inside if <= 1)',ylim=(.01,400))
axes[1].set(title='Conditional endpoint error-box coverage',ylabel='Observed error / certified local radius',ylim=(.001,30))
fig.suptitle('Second-round endpoint N10 trajectories: guarantee coverage gaps\nSaved-state checks; local boxes do not supply a priori closed-loop error bounds',fontsize=12)
fig.savefig(H/'guarantee_coverage.png',dpi=180);fig.savefig(H/'guarantee_coverage.svg');plt.close(fig)
print('FIGURES_COMPLETE')
