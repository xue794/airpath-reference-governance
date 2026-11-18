"""Independent saved-result validation, safe to rerun while cases finish.

No long waits or simulation launches. Missing/incomplete files are reported as
pending; scientific check failures are retained, not discarded.
"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / 'full6dof'
METHODS = ('T', 'N0', 'N10')
CASES = {
    'V31_C60': {'coarse': 'V31_C60_warm', 'fine': 'V31_C60_warm_fine'},
    'V32_C60': {'coarse': 'V32_C60_warm', 'fine': 'V32_C60_warm_fine'}
}
MARGIN_NAMES = ('reference_lower', 'reference_upper', 'feedback_lower', 'feedback_upper')
TOL = 1e-6


def same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    try:
        return bool(np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True))
    except TypeError:
        return bool(np.array_equal(np.asarray(a), np.asarray(b)))


def value_range(x):
    x = np.asarray(x)
    return dict(min=float(np.min(x)), max=float(np.max(x)))


def run_check(label, method, data):
    t = np.atleast_1d(data['Ts'])
    x, d, s = data['Xs'], data['D'], data['Summary']
    event_times = np.atleast_1d(s['event_times_s']).astype(float)
    event_indices = np.atleast_1d(s['event_indices']).astype(int)
    accepted_after_event = int(np.sum(t > event_times[0]+1e-8)) if len(event_times) else 0
    endpoint = s['status'] == 'terminal_original_root_boundary'
    completed = s['status'] == 'completed'
    horizon = float(s['requested_time_s'])
    check = dict(
        finite_state_and_diagnostics=bool(np.all(np.isfinite(x)) and np.all(np.isfinite(d))),
        strictly_increasing_time=bool(np.all(np.diff(t)>0)),
        all_original_root_margins_at_least_minus_1e_minus_6=bool(np.all(d[:,12:16] >= -TOL)),
        no_output_after_first_event=accepted_after_event == 0,
        completed_only_at_full_horizon=not completed or (abs(t[-1]-horizon)<1e-8 and len(event_times)==0),
        event_status_has_matching_last_output=not endpoint or (len(event_times)>0 and abs(t[-1]-event_times[0])<1e-8),
        all_conventional_state_checks_pass=bool(np.all(d[:,18:23]>0.5)),
        all_psi1_at_least_minus_1e_minus_6=bool(np.all(d[:,17]>=-TOL)),
        no_actual_thrust_box_violation=bool(np.all((d[:,11]>=-TOL)&(d[:,11]<=60+TOL)))
    )
    return dict(
        label=label, method=method, status=s['status'], samples=len(t),
        start_time_s=float(t[0]), end_time_s=float(t[-1]), requested_end_time_s=horizon,
        max_step_s=float(data['C']['MaxStep']), event_times_s=event_times.tolist(), event_indices=event_indices.tolist(),
        event_channels=[MARGIN_NAMES[i-1] for i in event_indices],
        accepted_points_after_first_event=accepted_after_event,
        minimum_root_margins_N=dict(zip(MARGIN_NAMES, np.min(d[:,12:16],axis=0).tolist())),
        negative_root_margin_samples=int(np.sum(np.any(d[:,12:16]<-TOL,axis=1))),
        minimum_psi1=float(np.min(d[:,17])),
        alpha_reference_deg=value_range(d[:,6]), alpha_command_deg=value_range(d[:,7]), alpha_actual_deg=value_range(d[:,8]),
        reference_thrust_N=value_range(d[:,9]), raw_command_thrust_N=value_range(d[:,10]), actual_thrust_N=value_range(d[:,11]),
        reference_thrust_reserve_violation_samples=int(np.sum((d[:,9]<3-TOL)|(d[:,9]>57+TOL))),
        raw_command_thrust_box_violation_samples=int(np.sum((d[:,10]<-TOL)|(d[:,10]>60+TOL))),
        maximum_QP_primal_normalized=float(np.max(d[:,23])),
        max_VB_over_rho=float(np.max(d[:,16])),
        extension_evaluations=s['extension_evaluations'], checks=check, all_checks_pass=all(check.values())
    )


def refinement(case, method, a, b):
    ta, tb = np.atleast_1d(a['Ts']), np.atleast_1d(b['Ts'])
    times, ia, ib = np.intersect1d(ta,tb,return_indices=True)
    difference = np.abs(a['Xs'][ia]-b['Xs'][ib])
    ca, cb = dict(a['C']), dict(b['C'])
    ca.pop('MaxStep'); cb.pop('MaxStep')
    sa, sb = a['Summary'], b['Summary']
    return dict(case=case,method=method,
        same_status=sa['status']==sb['status'], same_event_indices=same(sa['event_indices'],sb['event_indices']),
        end_time_difference_s=float(abs(ta[-1]-tb[-1])),
        coarse_end_time_s=float(ta[-1]),fine_end_time_s=float(tb[-1]),
        common_exact_saved_times=len(times),common_grid_last_time_s=float(times[-1]),
        max_state_abs_difference=float(np.max(difference)),
        max_state_abs_difference_by_component=np.max(difference,axis=0).tolist(),
        max_position_vector_difference_m=float(np.max(np.linalg.norm(a['Xs'][ia,:3]-b['Xs'][ib,:3],axis=1))),
        max_position_component_difference_m=float(np.max(difference[:,:3])),
        same_C_except_MaxStep=same(ca,cb),same_request_O=same(a['O'],b['O']),
        coarse_min_reference_T_N=float(np.min(a['D'][:,9])),fine_min_reference_T_N=float(np.min(b['D'][:,9])),
        coarse_min_feedback_upper_N=float(np.min(a['D'][:,15])),fine_min_feedback_upper_N=float(np.min(b['D'][:,15])))


def common_interval(case, level, data):
    start=max(float(np.atleast_1d(data[m]['Ts'])[0]) for m in METHODS)
    end=min(float(np.atleast_1d(data[m]['Ts'])[-1]) for m in METHODS)
    nodes=[np.array([start,end])]
    for m in METHODS:
        t=np.atleast_1d(data[m]['Ts']); nodes.append(t[(t>start)&(t<end)])
    grid=np.unique(np.concatenate(nodes))
    values=[]
    for m in METHODS:
        t=np.atleast_1d(data[m]['Ts']);e=data[m]['D'][:,1]
        interp=np.interp(grid,t,e)
        integral=float(np.trapezoid(interp**2,grid))
        values.append(dict(method=m, common_interval_RMS_m=float(np.sqrt(integral/(end-start))),
                           error_squared_integral_m2s=integral,
                           final_endpoint_interpolated=not bool(np.any(t==end)),
                           individual_end_time_s=float(t[-1]),individual_status=data[m]['Summary']['status']))
    return dict(case=case,level=level,start_time_s=start,end_time_s=end,duration_s=end-start,
                integration_nodes=len(grid),methods=values,
                calculation='Linearly interpolate saved position-error norm onto union grid plus exact common endpoint, then trapezoid(error^2)/duration. No extrapolation or zero padding.')


data={};pending=[];failures=[];checks=[]
for case,levels in CASES.items():
    data[case]={}
    for level,label in levels.items():
        data[case][level]={}
        for method in METHODS:
            file=ROOT/label/(method+'.mat')
            if not file.exists():
                pending.append(str(file));continue
            try:
                d=loadmat(file,simplify_cells=True)
                c=run_check(label,method,d)
                data[case][level][method]=d;checks.append(c)
                if not c['all_checks_pass']:
                    failures.append(dict(label=label,method=method,failed_checks=[k for k,v in c['checks'].items() if not v]))
            except Exception as e:
                pending.append(f'{file}: read incomplete/error: {e}')
refinements=[];intervals=[]
for case,levels in data.items():
    for method in METHODS:
        if all(method in levels[level] for level in ('coarse','fine')):
            r=refinement(case,method,levels['coarse'][method],levels['fine'][method]);refinements.append(r)
            if not r['same_status'] or not r['same_event_indices']:
                failures.append(dict(case=case,method=method,refinement_termination_changed=True))
    for level,runs in levels.items():
        if all(m in runs for m in METHODS):intervals.append(common_interval(case,level,runs))

report=dict(scope='Independent saved-state and common-time validation; numerical evidence, not a continuous-time safety proof.',
            all_expected_results_present=not pending, pending=pending,
            all_available_run_checks_pass=not failures,failures=failures,
            run_checks=checks,refinement=refinements,common_interval_metrics=intervals,
            limitations=[
                'No comparisons of per-method RMS over unequal stopping horizons.',
                'Root-margin event is the original selected inversion interval, not the full physical flight envelope.',
                'Sampled accepted states and MaxStep refinement do not prove between-sample invariance.',
                'State difference components have mixed units; use component-wise values, not an unscaled physical norm.'
            ])
(HERE/'final_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,
    default=lambda v: v.item() if isinstance(v,np.generic) else v.tolist()),encoding='utf-8')
lines=['最终结果独立核验（可重复执行）','',
       f"预期数据是否齐全：{not pending}；现有运行检查是否通过：{not failures}",
       '检查包括保存状态有限性、时间递增、四个原始根裕度、事件后无输出、完成时长、普通状态/psi1和实际推力盒。','']
for c in checks:
    lines.append(f"{c['label']} / {c['method']}: {c['status']}, t_end={c['end_time_s']:.12g}s, event={c['event_channels']}, samples={c['samples']}, pass={c['all_checks_pass']}")
    lines.append(f"  原始根最小裕度N={c['minimum_root_margins_N']}; reference/raw/actual推力N={c['reference_thrust_N']}/{c['raw_command_thrust_N']}/{c['actual_thrust_N']}")
    lines.append(f"  reference/command/actual迎角deg={c['alpha_reference_deg']}/{c['alpha_command_deg']}/{c['alpha_actual_deg']}")
lines+=['','粗/细MaxStep比较（同科学工况、同方法）：']
for r in refinements:
    lines.append(f"{r['case']} / {r['method']}: 同终止状态={r['same_status']}, 同事件={r['same_event_indices']}, 终止时间差={r['end_time_difference_s']:.12g}s;")
    lines.append(f"  精确共有采样点={r['common_exact_saved_times']}; 最大位置向量差={r['max_position_vector_difference_m']:.12g}m; 各状态分量差={r['max_state_abs_difference_by_component']}")
lines+=['','共同时间区间RMS（误差平方梯形积分除以共同持续时间再开方）：']
for c in intervals:
    lines.append(f"{c['case']} / {c['level']}: [{c['start_time_s']:.12g}, {c['end_time_s']:.12g}]s")
    for m in c['methods']:
        lines.append(f"  {m['method']}: RMS={m['common_interval_RMS_m']:.12g}m; 末端插值={m['final_endpoint_interpolated']}; 原轨迹终点={m['individual_end_time_s']:.12g}s")
lines+=['','不同方法终止时域不同，禁止用各自RMS_until_stop直接排名。这里只在同一个共同区间内积分；终止原因和是否完成35s另列。',
        '所有检查仍是数值证据；没有KB的方案不能继承原递归可行定理，0–9deg区间边界也不是飞机完整能力边界。']
if failures:lines+=['','检查失败或终止分类变化：',json.dumps(failures,ensure_ascii=False,indent=2)]
if pending:lines+=['','等待文件（本脚本不等待、不启动仿真）：']+pending
(HERE/'final_validation.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(dict(all_expected_results_present=not pending,available_runs=len(checks),pending_count=len(pending),
                     all_available_checks_pass=not failures,report=str(HERE/'final_validation.txt')),ensure_ascii=False))
