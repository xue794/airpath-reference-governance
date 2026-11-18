"""Assemble coverage, verified negative witnesses, and bounded search outcomes."""
from pathlib import Path
import json,sys
from fractions import Fraction
import numpy as np
from scipy.io import loadmat
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE/'original_certificate'))
import interval_certificate as ic
audit=json.loads((HERE/'backup_probe.json').read_text(encoding='utf-8'))
search=json.loads((HERE/'recentered_certificates.json').read_text(encoding='utf-8'))
follow=json.loads((HERE/'followup_certificates.json').read_text(encoding='utf-8'))
R=loadmat(HERE/'V31_backup_saved_audit.mat',simplify_cells=True)['R'];P=R['PB']
candidates=[dict(label='original_38.1_rho0.1',center=38.1,turn=0.,rho=.1,
    certificate_source='original_certificate/backup_rho_0p10_certificate.json')]
for r in search+follow:
    if r['certified']:
        center=float(r['center_mps']);turn=float(r.get('course_rate_center_rad_s',0));rho=float(r['rho'])
        candidates.append(dict(label=f'V{center:g}_turn{turn:g}_rho{rho:g}',center=center,turn=turn,rho=rho))
best={}
for c in candidates:
    key=(c['center'],c['turn'])
    if key not in best or c['rho']>best[key]['rho']:best[key]=c
candidates=list(best.values())

def distance(x,c):
    z=x[:,[0,1,3,4,5]].copy();z[:,0]-=c['center'];z[:,4]-=c['turn']
    return np.einsum('ni,ij,nj->n',z,P,z)

datasets={}
for v in [31,32]:
    d=loadmat(HERE/f'V{v}_backup_saved_audit.mat',simplify_cells=True)
    datasets[f'full6dof_V{v}_N10']=(np.atleast_1d(d['t']),d['X'])
r=np.genfromtxt(ROOT/'endpoint_probe_20261106/reference/continuous_full17_N10_maxstep0.005.csv',delimiter=',',names=True)
datasets['reference_only_V31_N10']=(r['t'],np.column_stack([r[k] for k in ('V','gamma','chi','nuV','nuGamma','nuChi')]))
coverage=[]
for name,(t,x) in datasets.items():
    ratios=[];records=[]
    for c in candidates:
        vb=distance(x,c);ratio=vb/c['rho'];ratios.append(ratio)
        outside=np.flatnonzero(ratio>1+1e-10)
        records.append(dict(**c,covered_saved_points=int(np.sum(ratio<=1+1e-10)),samples=len(t),
            minimum_VB=float(np.min(vb)),required_rho_for_saved_trajectory=float(np.max(vb)),
            max_required_over_certified_rho=float(np.max(ratio)),
            first_outside_time_s=float(t[outside[0]]) if len(outside) else None))
    ratios=np.column_stack(ratios);union=np.any(ratios<=1+1e-10,axis=1)
    keypoints=[]
    for target in [0,10,14.65,15.4,16.8,17.52,19,25,35]:
        k=int(np.argmin(abs(t-target)));j=int(np.argmin(ratios[k]))
        keypoints.append(dict(requested_time_s=target,saved_time_s=float(t[k]),x=x[k].tolist(),
            best_candidate=candidates[j]['label'],best_VB_over_rho=float(ratios[k,j]),
            VB_over_rho_by_candidate={c['label']:float(ratios[k,jj]) for jj,c in enumerate(candidates)}))
    coverage.append(dict(dataset=name,samples=len(t),candidates=records,
        union_covered_saved_points=int(np.sum(union)),
        union_coverage_is_not_a_switching_certificate=True,keypoints=keypoints))

# Verify counterexample signs with exact rational coordinates from saved doubles.
# These concern the explicitly transcribed nominal reference model, not an
# interval enclosure of the numerical ODE solution as a function of time.
row_keys={'N_plus':'phiNplus','T_min':'phiTmin','T_max':'phiTmax'}
witnesses=[]
for c in audit['cases']:
    original=[v for v in c['centers'] if v['Vb']==38.1][0]
    for row,key in row_keys.items():
        j=audit['row_names'].index(row);xf=original['row_min_x'][j]
        v,gam,chi,nv,ng,nc=xf;ic.Vb=ic.I('38.1')
        exact=lambda f:ic.I(str(Fraction.from_float(float(f))))
        z=[exact(v)-ic.Vb,exact(gam),exact(nv),exact(ng),exact(nc)]
        checks,valid=ic.evaluate(z);vb=ic.ellipse(z)
        assert valid and checks[key].hi<0
        witnesses.append(dict(case=c['speed_request'],row=row,time_s=original['row_min_time_s'][j],
            saved_x=xf,backup_u=original['row_min_u'][j],MATLAB_value=original['row_min'][j],
            interval_value=checks[key].pair(),strict_negative_verified=True,
            Tff_interval=(checks['hTmin']+3).pair(),hNplus_interval=checks['hNplus'].pair(),
            VB_interval=vb.pair(),root_and_10N_reference_reserve_positive=valid))

static=np.array([[31.,0,0,0,0,0]])
result=dict(scope='Nominal reference backup research. Exact interval certificates are separated from saved-point coverage and numerical trajectory evidence.',
    original_backup_audit=audit,largest_certified_candidates_in_this_finite_search=candidates,
    candidates_are_not_maximal=True,coverage=coverage,verified_negative_witnesses=witnesses,
    disclosed_static_endpoint_witness_VB={c['label']:float(distance(static,c)[0]) for c in candidates},
    bounded_search_results=search+follow,
    local_obstructions=json.loads((HERE/'local_obstructions.json').read_text(encoding='utf-8')),
    kernel_regression=json.loads((HERE/'kernel_regression.json').read_text(encoding='utf-8')),
    limitations=['Fixed backup failure does not imply nonexistence of another feasible backup.',
        'A union of local invariant sets is not a certified recovery/switching construction without transition checks.',
        'No finite-horizon recovery or trajectory tube has been certified in this bounded study.',
        'The original interval kernel certifies the explicitly transcribed nominal reference model; it is not a formal proof of MATLAB execution or aircraft safety.'])
(HERE/'backup_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')

lines=['第三轮：备份保证与第二轮收益区域的覆盖核查','',
 '结论：原方案不能直接把第二轮无 KB 的收益轨迹纳入原递归可行定理。',
 '原因不仅是椭球覆盖不足；原固定备份在若干收益状态上也明确违反操作 CBF 行。',
 '但这不否定其他备份构造存在。已用原外向区间核验证了新的非零局部邻域。','',
 '一、严格局部参考证书（本有限搜索中的最大已通过值，并非全局最大）']
for c in candidates:
    lines.append(f"Vb={c['center']:g} m/s, nuChiB={c['turn']:g} rad/s, rho={c['rho']:g}")
lines += ['', '原 38.1/.1 证书沿用已归档证据；速度中心候选复用原整数外向算术与原 evaluate。',
 '转弯候选使用明确修改后的核，源码生成器逐处断言且完整保存，并非宣称原核完全未改：',
 '修改物理 course rate 与误差 rate 的关系、备份 uChi，以及两条 course-rate CBF 常数。',
 '平移速度中心不改变 A_B、Q_B、P_B。恒转弯候选以 nuChi-nuChiB 为第五误差坐标，',
 '并用 uChi=-1.5*(nuChi-nuChiB)，其余反馈与 CBF 增益不变；这是 course 无关',
 '参考坐标上的相对平衡点，不是实际飞机稳态转弯证书。备份行解析正下界为 .145*rho。',
 '所有局部证书含原 10 N 端点预留、3–57 N 推力预留、初始化、正则性及 1:17 行。',
 '转弯核在零转率下通过 20 个与原核逐整数端点一致的检查；在非零转率下',
 '通过 11 个原 MATLAB 仿射核对照点，最大区间外数值差约 3.56e-15。',
 '这些有限回归只核对转译，不替代另存的全椭球区间覆盖证书。','',
 '二、原备份在完整 6DOF 轨迹所保存参考状态上的结果']
for case in audit['cases']:
    c=[c for c in case['centers'] if c['Vb']==38.1][0]
    lines.append(f"{case['speed_request']} m/s: 原 KB 覆盖 {c['covered_samples_rho01']}/{case['samples']} 点；max VB={c['maxVB']:.9f}，相对原 rho=.1 为 {c['maxVB']/.1:.3f} 倍；原备份满足全部 1:17 行 {c['backup_all17_pass_samples']}/{case['samples']} 点。")
for w in witnesses:
    lines.append(f"{w['case']} m/s, t={w['time_s']:.2f}s, 原备份 {w['row']} 行={w['MATLAB_value']:.9f} N/s，外向区间负号已独立确认；Tff∈{w['Tff_interval']} N。")
lines += ['', '以上违反的是在这些状态上假设采用原备份后的导数条件；保存轨迹实际使用的是第二轮 N10 QP，',
 '不能误写为实际执行了不合法备份，也不能由此推出飞机物理失效。',
 '把中心改为 31、32、34 m/s 但保持一个固定备份，也没有让该备份在全部保存点可行；完整逐行记录见 JSON。','',
 '三、严格局部集合与轨迹覆盖']
for ds in coverage:
    lines.append(f"{ds['dataset']}: 已认证候选集合的并集覆盖保存点 {ds['union_covered_saved_points']}/{ds['samples']}；并集数字不构成切换/恢复证书。")
    for c in ds['candidates']:
        lines.append(f"  {c['label']}: 覆盖 {c['covered_saved_points']}/{c['samples']}，所需 max VB={c['required_rho_for_saved_trajectory']:.6f}，为本候选已认证 rho 的 {c['max_required_over_certified_rho']:.1f} 倍。")
    for k in ds['keypoints']:
        if k['requested_time_s'] in (14.65,17.52,19,25,35):
            lines.append(f"  t={k['saved_time_s']:.2f}s，距离最近的已认证候选 {k['best_candidate']}，VB/rho={k['best_VB_over_rho']:.3f}。")
lines += ['', '31 m/s 水平、零参考率的已披露静态端点例，现在位于 31 m/s 局部证书中心。',
 '这支持在收益相关低速区域建立非零局部证书的可行性，但不能把该中心的静态机制证据',
 '替代完整减速—转弯—恢复轨迹覆盖，也不能据此断言端点行在包含于根域的 KB 内有独立状态保护收益。','',
 '四、未通过候选与结论边界']
for r in search+follow:
    if not r['certified']:
        lines.append(f"Vb={r['center_mps']}, rho={r['rho']}: 限定节点数内未完成区间覆盖，nodes={r['nodes']}，最后未解决行={list(r.get('failed',{}))}；单凭该状态不能称为不可行。")
for r in result['local_obstructions']:
    if r['negative_sign_and_interior_strictly_verified']:
        lines.append(f"另有严格内部反例：Vb={r['center']}, rho={r['rho']}，{r['target']}∈{r['target_interval']}；该点否定整个候选椭球满足所需预留，而不是优化器判据。")
lines += ['原固定形状/增益的简单扩大已经有明确阻碍；未搜索或否定所有其他形状和反馈。',
 '本轮未认证轨迹备份管、有限时域恢复或候选之间的转移，因而不能把局部集合并集接成新的递归可行主张。',
 '下一步若保留第二轮动态收益，需要设计与状态/任务相适应的恢复备份并认证过渡区域，',
 '或把论文主张限定为能被局部证书覆盖的任务，并重新验证该任务的实际收益。','',
 '复现：MATLAB 执行 run_backup_probe、run_turn_crosscheck；Python 依次执行',
 'certify_recentered.py、certify_turn_center.py、search_local_obstruction.py、',
 'verify_local_kernel.py、summarize_backup.py。',
 '原两轮数据与代码只读；本目录包含所有未通过候选和中间检查记录。']
(HERE/'BACKUP_RESULTS_CN.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(dict(summary=str(HERE/'backup_summary.json'),candidates=len(candidates),verified_negative_witnesses=len(witnesses)),ensure_ascii=False))
