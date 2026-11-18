"""Independent comparison of cold and warm numerical QP starts.

The compared mathematical QP is unchanged. This checks saved full states,
events and diagnostics, not just scalar tracking summaries.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / 'full6dof'


def same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    try:
        return bool(np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True))
    except TypeError:
        return bool(np.array_equal(np.asarray(a), np.asarray(b)))


snapshot = (ROOT / 'ep6_snapshot.m').read_text()
expected = snapshot.replace('function S=ep6_snapshot(', 'function S=ep6_snapshot_warm(', 1)
expected = expected.replace('Gov=ep6_governor(', 'Gov=ep6_governor_warm(', 1)
assert expected == (ROOT / 'ep6_snapshot_warm.m').read_text()

records = []
for method in ('T', 'N0', 'N10'):
    oldfile = ROOT / 'V31_C60' / (method + '.mat')
    newfile = ROOT / 'V31_C60_warm' / (method + '.mat')
    a = loadmat(oldfile, simplify_cells=True)
    b = loadmat(newfile, simplify_cells=True)
    ta, tb = np.atleast_1d(a['Ts']), np.atleast_1d(b['Ts'])
    common, ia, ib = np.intersect1d(ta, tb, return_indices=True)
    delta = np.abs(a['Xs'][ia] - b['Xs'][ib])
    sa, sb = a['Summary'], b['Summary']
    by_component = np.max(delta, axis=0)
    rec = dict(
        method=method,
        original_mat_sha256=hashlib.sha256(oldfile.read_bytes()).hexdigest(),
        warm_mat_sha256=hashlib.sha256(newfile.read_bytes()).hexdigest(),
        original_status=sa['status'], warm_status=sb['status'],
        parameter_structure_exact_match=same(a['C'], b['C']),
        request_options_exact_match=same(a['O'], b['O']),
        initial_state_record_exact_match=same(a['I0'], b['I0']),
        original_samples=len(ta), warm_samples=len(tb),
        common_exact_times=len(common),
        end_time_difference_s=float(abs(ta[-1] - tb[-1])),
        original_event_indices=np.atleast_1d(sa['event_indices']).tolist(),
        warm_event_indices=np.atleast_1d(sb['event_indices']).tolist(),
        max_state_abs_difference=float(np.max(delta)),
        max_difference_by_state_component=by_component.tolist(),
        max_position_component_difference_m=float(np.max(delta[:, :3])),
        max_position_vector_difference_m=float(np.max(np.linalg.norm(a['Xs'][ia, :3] - b['Xs'][ib, :3], axis=1))),
        last_output_state_difference=np.abs(a['Xs'][-1]-b['Xs'][-1]).tolist(),
        absolute_RMS_difference_m=float(abs(sa['RMS_until_stop_m']-sb['RMS_until_stop_m'])),
        original_minimum_root_margin_by_channel_N=np.min(a['D'][:, 12:16], axis=0).tolist(),
        warm_minimum_root_margin_by_channel_N=np.min(b['D'][:, 12:16], axis=0).tolist(),
        original_saved_root_domain_violations=int(np.sum(np.any(a['D'][:, 12:16] < -1e-6, axis=1))),
        warm_saved_root_domain_violations=int(np.sum(np.any(b['D'][:, 12:16] < -1e-6, axis=1))),
        original_max_QP_primal_normalized=float(sa['max_QP_primal_normalized']),
        warm_max_QP_primal_normalized=float(sb['max_QP_primal_normalized']),
        original_runtime_s=float(sa['wall_seconds']), warm_runtime_s=float(sb['wall_seconds'])
    )
    assert rec['parameter_structure_exact_match'] and rec['request_options_exact_match'] and rec['initial_state_record_exact_match']
    assert rec['original_status'] == rec['warm_status']
    assert rec['original_event_indices'] == rec['warm_event_indices']
    records.append(rec)

result = dict(
    scope='Numerical-start equivalence check for V31_C60, MaxStep=0.01; not a proof for all future cases.',
    static_review='Only feasible-start proposal changed; original w2b_qp3 verifies current feasibility and KKT, retaining fallback.',
    snapshot_only_dispatch_changed=True,
    formula='As=A_W diag(scale); z=As^T(As As^T)^(-1)(b_W-A_W u0); u_start=u0+diag(scale)z',
    tolerance_interpretation='T/N0 agree at near-roundoff scale; N10 differs by micrometers in position, with unchanged termination and qualitative conclusions.',
    records=records)
(HERE/'warm_start_audit.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')

lines = [
    'Warm-start 独立实现和轨迹审计（2026-11-06）', '',
    '结论：改变的是数值可行起点提议，未改变数学QP。原求解器继续验证当前约束和KKT，保留枚举回退。',
    '当前V31_C60粗步长对照数值一致，但N10不是逐位相同。可以使用warm版本继续同问题步长细化，',
    '停止完全冗余的原版细化计算；不得把不同科学工况混为冗余。', '',
    '静态代码核查：',
    '- A、b、Scale、端点/推力行集合、DeltaN设定以及w2b_qp3均未变化。',
    '- 上次正对偶支撑行仅产生新的加权等式投影起点；公式正确。',
    '- 只有同method复用支撑；旧ODE试探步的支撑也要接受当前全部约束检查。',
    '- 若提议不满足约束，仍用其他起点或原枚举流程；没有把历史控制直接当成解。',
    '- ep6_snapshot_warm与ep6_snapshot的差别逐字检查，仅函数名和governor调用名变化。', '',
    '逐轨迹数值核验：'
]
for r in records:
    lines += [
        f"{r['method']}: 状态={r['warm_status']}; 原/新样本={r['original_samples']}/{r['warm_samples']}; "
        f"共同精确时间点={r['common_exact_times']}; 终止时间差={r['end_time_difference_s']:.12g}s;",
        f"  全36维最大绝对差={r['max_state_abs_difference']:.12g}（混合单位，仅作逐项数值一致性）；"
        f"位置分量最大差={r['max_position_component_difference_m']:.12g}m; RMS差={r['absolute_RMS_difference_m']:.12g}m;",
        f"  原/新最大归一化primal残差={r['original_max_QP_primal_normalized']:.6g}/{r['warm_max_QP_primal_normalized']:.6g}; "
        f"保存状态原始根域违例={r['original_saved_root_domain_violations']}/{r['warm_saved_root_domain_violations']};",
        f"  C、O、I0参数记录完全相等；原/新时间成本={r['original_runtime_s']:.3f}/{r['warm_runtime_s']:.3f}s（并行环境下非正式计时基准）。"
    ]
lines += ['',
    '解释边界：',
    '数学目标严格凸，改变起点不改变唯一最优解。浮点求解路径改变可导致细小轨迹差异，',
    '不能声称warm与原实现逐位完全等价。当前N10位置差为微米量级，远小于报告的米级误差',
    '和数牛顿端点裕度，不改变科学结论。后续新工况仍应报告求解残差和域检查。',
    '不同方法若终止时域不同，其RMS不能直接比较；本审计只比较同method原/新实现。'
]
(HERE/'warm_start_audit.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
print(json.dumps(dict(all_complete=True, records=len(records), report=str(HERE/'warm_start_audit.txt'))))
