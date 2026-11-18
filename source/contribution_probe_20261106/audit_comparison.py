"""Independent saved-trajectory and refinement audit for the AoA comparisons."""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat

H=Path(__file__).resolve().parent
F=H/'full6dof'
OLD=H.parent/'endpoint_probe_20261106'/'full6dof'
expected=[(31,'A0',.01),(31,'Am',.01),(31,'A0',.005),(31,'Am',.005),
          (32,'Am',.01),(32,'Am',.005)]
runs={}; checks=[]; pending=[]
for speed,method,step in expected:
    tag=f'V{speed}_C60'+('_fine' if step==.005 else '')
    p=F/tag/f'{method}.mat'
    if not p.exists():pending.append(str(p));continue
    m=loadmat(p,simplify_cells=True)
    t=np.asarray(m['Ts']).reshape(-1);x=np.asarray(m['Xs']);d=np.asarray(m['D'])
    s=m['Summary'];te=np.asarray(s['event_times_s']).reshape(-1)
    tests={'finite_states':bool(np.isfinite(x).all()),'increasing_time':bool((np.diff(t)>0).all()),
           'finite_diagnostics':bool(np.isfinite(d).all()),
           'root_margins':bool((d[:,12:16]>=-1e-6).all()),
           'conventional_reference':bool((d[:,18:23]==1).all()),
           'hocbf_initialization':bool(d[:,17].min()>=-1e-7),
           'reference_thrust':bool(d[:,9].min()>=3-1e-6 and d[:,9].max()<=57+1e-6),
           'hardware_thrust':bool(d[:,10:12].min()>=-1e-6 and d[:,10:12].max()<=60+1e-6),
           'qp_primal':bool(d[:,23].max()<1e-8),
           'completed_or_original_event':s['status'] in ['completed','terminal_original_root_boundary']}
    if method=='Am':tests['matched_10N_reference_set']=bool(d[:,12:14].min()>=10-1e-6)
    if s['status']=='completed':tests['completed_duration']=bool(abs(t[-1]-35)<1e-9)
    if len(te):tests['no_post_event_output']=bool(t[-1]<=te[0]+1e-8 and abs(t[-1]-te[0])<1e-8)
    r={'case':tag,'method':method,'step':step,'passed':all(tests.values()),'checks':tests,
       'status':s['status'],'end_time_s':float(t[-1]),'min_reference_margin_N':float(d[:,12:14].min()),
       'min_feedback_margin_N':float(d[:,14:16].min()),'min_reference_thrust_N':float(d[:,9].min()),
       'max_reference_thrust_N':float(d[:,9].max()),'max_QP_primal':float(d[:,23].max()),
       'position_RMS_time_integral_until_stop_m':float(np.sqrt(np.trapezoid(d[:,1]**2,t)/(t[-1]-t[0])))}
    checks.append(r);runs[speed,method,step]=(t,x,d,s)
refinements=[]
for speed,method in [(31,'A0'),(31,'Am'),(32,'Am')]:
    if (speed,method,.01) not in runs or (speed,method,.005) not in runs:continue
    a=runs[speed,method,.01];b=runs[speed,method,.005]
    common,ia,ib=np.intersect1d(np.round(a[0],10),np.round(b[0],10),return_indices=True)
    dif=np.abs(a[1][ia]-b[1][ib])
    refinements.append({'speed':speed,'method':method,'same_status':a[3]['status']==b[3]['status'],
                        'same_events':bool(np.array_equal(a[3]['event_indices'],b[3]['event_indices'])),
                        'end_time_difference_s':float(abs(a[0][-1]-b[0][-1])),
                        'common_saved_points':len(common),'max_state_component_difference':dif.max(axis=0).tolist(),
                        'max_position_vector_difference_m':float(np.linalg.norm(dif[:,:3],axis=1).max())})
pairs=[]
for speed,am,endpoint in [(31,'A0','N0'),(31,'Am','N10'),(32,'Am','N10')]:
    if (speed,am,.005) not in runs:continue
    aa=runs[speed,am,.005];em=loadmat(OLD/f'V{speed}_C60_warm_fine/{endpoint}.mat',simplify_cells=True)
    et=np.asarray(em['Ts']).reshape(-1);ed=np.asarray(em['D']);end=min(et[-1],aa[0][-1])
    grid=np.unique(np.r_[et[et<=end],aa[0][aa[0]<=end],end])
    q={}
    for name,t,d in [(am,aa[0],aa[2]),(endpoint,et,ed)]:
        err=np.interp(grid,t,d[:,1]);q[name]={'common_interval_RMS_m':float(np.sqrt(np.trapezoid(err**2,grid)/(grid[-1]-grid[0]))),
          'individual_end_time_s':float(t[-1]),'max_position_error_until_stop_m':float(d[:,1].max()),
          'min_reference_upper_N':float(d[:,13].min()),'min_feedback_upper_N':float(d[:,15].min()),
          'max_VB_over_rho':float(d[:,16].max())}
    pairs.append({'speed':speed,'methods':[am,endpoint],'common_end_s':float(end),'metrics':q})
result={'complete':not pending,'pending':pending,'all_available_run_checks_pass':all(c['passed'] for c in checks),
        'run_checks':checks,'refinement':refinements,'matched_comparisons':pairs,
        'limitations':['Saved-point checks are not interval proofs.','No KB imposed; original recursive-feasibility theorem does not transfer.',
                       'AoA and endpoint sets are matched, fixed CBF gains need not give identical interior derivative rows.']}
(H/'comparison_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({'complete':result['complete'],'available':len(checks),'pending':len(pending),
                  'checks_pass':result['all_available_run_checks_pass'],'comparisons':pairs},ensure_ascii=False))
