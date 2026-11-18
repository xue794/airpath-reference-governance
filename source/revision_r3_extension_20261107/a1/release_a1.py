"""Build the readable A1 handoff and an integrity-checked, self-contained release."""
import csv
import hashlib
from html import escape
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def read(name):
    with (HERE/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def html_table(head,rows):
    title=''.join('<th>'+escape(x)+'</th>' for x in head)
    body=''.join('<tr>'+''.join('<td>'+escape(str(c))+'</td>' for c in r)+'</tr>' for r in rows)
    return '<div class="scroll"><table><thead><tr>'+title+'</tr></thead><tbody>'+body+'</tbody></table></div>'
def fnum(r,k,n=4):return f'{float(r[k]):.{n}f}'
def span(data):return ' 至 '.join(f'{v:.4f}' for v in data)

def report(z,runs,ev):
    case_rows=[(r['id'],r['full_old_m'],r['full_abs_diff_m'],r['late_abs_diff_m'],r['support_times_s'],
                fnum(r,'attitude_peak_deg'),r['all_18_support_binding_match']) for r in runs]
    event_rows=[(r['id'],fnum(r,'time_s',2),fnum(r,'upper_dual_scaled',6),r['backup_support'],
                 fnum(r,'error_box_mps'),fnum(r,'attitude_sum_deg'),fnum(r,'width_J_mps'),
                 r['E_included']+'/'+r['G_included']+'/'+r['J_included'],fnum(r,'fb_upper_N'),fnum(r,'Traw_N')) for r in ev]
    ids_out=[r['id']+' @ '+fnum(r,'time_s',2)+' s' for r in ev if not int(r['J_included'])]
    joint_details='，'.join(ids_out) if ids_out else '全部入盒'
    h="""<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A1 原始轨迹复现与端点事件核查</title><style>
body{font:16px/1.75 "Microsoft YaHei",sans-serif;color:#172830;max-width:1060px;margin:auto;padding:30px min(4vw,50px) 65px}h1{font-size:30px;line-height:1.25}h2{font-size:21px;color:#185165;margin:32px 0 8px}p{margin:8px 0}
table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}td,th{text-align:right;border-bottom:1px solid #cfdce1;padding:5px 8px;white-space:nowrap}th{background:#e7f1f2}td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}a{color:#145367}code{font-size:.92em}small{color:#45606b}
@media print{body{padding:0}thead{display:table-header-group}tr{break-inside:avoid}.scroll{overflow:visible}}
</style><main>"""
    h+=f'''<h1>A1 原始轨迹复现与端点事件核查</h1>
<p><strong>{z["accepted"]}/17 条 Matched-QP 轨迹通过预定复现门槛。</strong>保存网格共 {z["grid_points"]:,} 点。原上端点行的 {z["upper_support_samples"]} 个支持点全部恢复。筛选分母为原 431 条完成的 Matched-QP 轨迹，69 条退出记录不进入此回查。这是对已有角色台账的定向复现。</p>
<h2>整条轨迹的复现</h2>
<p>在 MATLAB R2018a 中核对了 66 份原始依赖、39 个发行包文件和 500 份冻结扰动记录。种子 20260819，内存中 ρB=0.10，原 ODE45 配置、公称 120 s、0.05 s 保存网格。每例从原初始化重新积分，保存 36 维状态、原始控制输出及逐点 18 行系数。</p>
<p>所有完成与 CorePass 标志一致。306 个“试验×约束行”逐一匹配支持计数和绑定计数。全程位置 RMS 与原双精度 MAT 记录的最大绝对差为 {z["rms_full_max_abs_diff_m"]:.6g} m，末 60 秒为 {z["rms_late_max_abs_diff_m"]:.6g} m。成功判据为每项差值不超过 max(10⁻⁶ m，原 RMS×10⁻⁸)。Python 另从 36 维状态、解析位置指令、逐例风场和留存的 QP 系数重建 RMS、航迹、约束余量与对偶角色。完整原试验仅留有汇总，因此复现依据为归档的逐例指标与计数；完整时序是本次新保存的数据。</p>
<p>最差归一化 KKT 残差 {max(float(r['kkt_max']) for r in runs):.4g}（门槛 5×10⁻⁸），最差归一化行违反量 {max(float(r['norm_row_violation_max']) for r in runs):.4g}（门槛 10⁻⁸），最大 V<sub>B</sub>/ρ<sub>B</sub> {max(float(r['VB_over_rho_max']) for r in runs):.12g}（门槛 1+10⁻⁸）。数值详细值见 CSV。</p>'''
    h+=html_table(['样本','原全程 RMS (m)','全程差 (m)','后段差 (m)','支持时间 (s)','姿态和峰值 (°)','18 行计数匹配'],case_rows)
    h+=f'''<h2>21 个原上端点支持时刻</h2>
<p>时刻范围 {z['min_event_time_s']:.2f}–{z['max_event_time_s']:.2f} s，均为保存网格。备份行同时具有正对偶支持 {z['coactive_backup']}/21。使用各自真实扰动样本重建航迹误差后，端点充分盒 E 包含 {z['E_inside']}/21，几何充分盒 G 包含 {z['G_inside']}/21。联合充分盒 J 的航迹误差条件满足 {z['J_air_inside']}/21，姿态误差和不超过 0.70° 为 {z['attitude_inside']}/21；两项共同满足 {z['J_inside']}/21。J 未包含的样本与时刻：{joint_details}。</p>
<p>观测等效航迹误差范围 {span(z['event_error_range_mps'])} m/s，姿态误差和 {span(z['event_attitude_range_deg'])}°。反馈侧原始上端点力裕度 {span(z['event_fb_upper_range_N'])} N，未裁剪反馈推力 {span(z['event_raw_thrust_range_N'])} N。两个原始反馈端点力均严格为正 {z['fb_root_bracket_observed']}/21，未裁剪推力严格处于 (0,60) N 共 {z['raw_thrust_inside']}/21。参考推力预留为 3–57 N，模型实际执行器推力另记。</p>
<p>误差宽度以 max(|eV|,38.1|eγ|,38.1|eχ|) 计，角以 rad 计。J 半宽沿用 WP3 的外向区间方法，每个 MC 状态重新计算。表中的 1/0 表示该保存点实际误差是否落入对应充分盒。盒外点仍可直接检查实际命令。冻结的 1 m/s 搜索上限适用于 E/G/J 的半宽，达到上限只证明已验证到该值。</p>'''
    exclusions=z['J_exclusions']
    h+=f'''<p>未进入 J 的原因按互斥条件计数：仅航迹误差盒条件未通过 {exclusions['air_only']} 点，仅姿态条件未通过 {exclusions['attitude_only']} 点，两项均未通过 {exclusions['both']} 点。原 QP 上端点<strong>导数行</strong>的最大绝对剩余量为 {z['max_abs_upper_phi_Nps']:.5g} N/s；相同时刻<strong>原始参考端点力间隙</strong>仍为 {span(z['event_ref_upper_range_N'])} N。导数行参与和当前根接近 9° 是不同的判断。</p>'''
    for r in ev:
        if not int(r['raw_thrust_inside']):
            h+=f'''<p><strong>样本 {r['id']}，{float(r['time_s']):.2f} s。</strong>参考推力 {float(r['Tff_N']):.4f} N，参考与反馈攻角分别为 {float(r['alpha_ref_deg']):.4f}°、{float(r['alpha_fb_deg']):.4f}°，两者均在 0–9° 反演区间内。此时反馈推力原始请求为 {float(r['Traw_N']):.4f} N，模型执行器裁剪目标为 {float(r['Ttarget_clipped_N']):.4f} N，实际推力状态为 {float(r['Tact_N']):.5g} N。该点在 J 盒外。原试验 CorePass 包含对模型实际执行器的检查。仅凭本次原轨迹不能把这次原始指令越界归因于端点约束。</p>'''
    group_rows=[]
    for name,key in [('初始时刻 t=0','initial_events'),('其余保存时刻 t>0','later_events')]:
        g=z[key]; group_rows.append((name,len(g['trial_ids']),g['samples'],g['backup_support'],g['E_inside'],g['G_inside'],g['attitude_inside'],g['J_inside']))
    h+=html_table(['按时刻分组','试验数¹','支持点','备份同时支持','E 包含','G 包含','姿态和达标','J 包含'],group_rows)
    h+='<p><small>¹ 同一试验可在两组出现。按保存时刻计数不能换算成约束连续作用的时长。</small></p>'
    h+=html_table(['样本','t (s)','归一化对偶','备份支持','误差 (m/s)','姿态和 (°)','J 半宽 (m/s)','E/G/J 包含','反馈上裕度 (N)','原始推力 (N)'],event_rows)
    h+='''<h2>交接与数据口径</h2><p>上端点行的正对偶支持、与备份行共同支持、盒内反馈裕度均可由 A1 核实。端点约束是否独立改变指令及后续闭环，需要下一阶段在同状态、同请求下删除对应约束再比较。上述分母限于指定试验与保存时刻。</p>
<p><a href="trial_reproduction.csv">17 条复现明细</a>　<a href="row_counts.csv">306 行支持与绑定对照</a>　<a href="event_checks.csv">21 点全部事件字段</a>　<a href="a1_verified.json">汇总 JSON</a>　<a href="README_复算说明.txt">复算与目录说明</a></p></main></html>'''
    (HERE/'A1_复现与事件核查.html').write_text(h,encoding='utf-8')

GUIDE="""Aerospace 大修 A1 复现与事件核查

打开 A1_复现与事件核查.html 可查看结论、全部 17 例和 21 个保存时刻。三份 CSV 以 UTF-8 BOM 编码，可直接打开。

目录与信号
runs/id_XXX.mat：原入口 w2e_simulate_trial 重跑的完整 S 结构，含 t、36 维 x、原 Y、C、初态和逐例运行时间，另存 sample、本次 E 和原检查点 old。存档的 old 为当年汇总。
runs/id_XXX_rows.mat：与 S.t 同序的 affineC [n,18]、affineG [n,18,3]、物理尺度 dual 与归一化 dualScaled、VBoverRho、参考/反馈攻角（rad）及按下/上排序的原始端点力间隙（N）。计算条件使用 Y.u0 与 Y.uGov。绑定依据原缩放行容差，正对偶支持依据冻结的归一化对偶门槛。
checks：原重跑及代数快照检查。identity：源路径、冻结输入、选样、worker 身份及哈希。logs 与 console.log：运行记录（原 MATLAB 部分系统消息使用 Windows 本地编码）。
trial_reproduction.csv 给出各试验双精度原值、新值、门槛量和支持时刻。row_counts.csv 包含逐试验逐行的原与新支持及绑定计数。event_checks.csv 区别原始端点间隙、导数行余量、未裁剪反馈推力及实际推力，分别为 N、N/s、N、N。width_E/G/J 为逐点充分半宽，单位 m/s。J_air_included 为单独航迹盒条件，J_included 同时要求态势误差和 ≤0.70°。

复核
1. 保持解压后的目录结构。MATLAB R2018a + Parallel Computing Toolbox 可执行并行原试验。Python 3.11 或更新版本，需 NumPy 与 SciPy，mpmath 1.3.0 已随 WP3 区间核提供。
2. 在 revision_r3_extension_20261107/a1 运行 python verify_a1.py。它读取各试验双精度轨迹，独立重建检查并生成 CSV/JSON。python release_a1.py 生成阅读报告及归档。
3. 如需重新做积分，请在单独工作副本按相同相对目录放入 run_a1.m 及冻结输入，让新副本的 checks、runs 等结果目录从空开始。在 MATLAB 切换至新副本的 a1，运行 run_a1。已存在的逐例结果标识会跳过该例，包括失败记录；如有失败，应保留首轮记录并在独立副本诊断。原 500 对主程序会启动整批，请使用本包的定向入口。核验结论适用于这 17 条 QP 轨迹，未执行删行试验。
"""

def main():
    z=json.loads((HERE/'a1_verified.json').read_text())
    runs=read('trial_reproduction.csv'); events=read('event_checks.csv')
    report(z,runs,events)
    (HERE/'README_复算说明.txt').write_text(GUIDE.replace('态势误差','姿态误差'),encoding='utf-8-sig')
    paths={p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.tmp','.pyc'} and p.name!='release_manifest.json'}
    records=json.loads((HERE/'identity/input_hashes.json').read_bytes())
    paths.update(ROOT/p for p in records)
    core=ROOT/'FWUAV/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel'
    paths.update(core.glob('*.m'))
    vendor=ROOT/'revision_wp3_20261107/vendor'
    paths.update(p for p in vendor.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(paths)}
    mfile=HERE/'release_manifest.json';mfile.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    paths.add(mfile)
    out=ROOT/'Aerospace_A1复现与事件核查材料包_v1_0.zip';partial=out.with_suffix('.partial.zip')
    with ZipFile(partial,'w',compression=ZIP_DEFLATED,compresslevel=3,allowZip64=True) as zf:
        for p in sorted(paths):zf.write(p,p.relative_to(ROOT).as_posix())
    with ZipFile(partial) as zf:
        assert set(zf.namelist())==set(manifest)|{mfile.relative_to(ROOT).as_posix()}
        for name,meta in manifest.items():
            with zf.open(name) as f:assert hashlib.file_digest(f,'sha256').hexdigest()==meta['sha256'],name
    partial.replace(out)
    print(len(paths),'verified entries;',out.stat().st_size,'bytes; SHA-256:',digest(out))
if __name__=='__main__': main()
