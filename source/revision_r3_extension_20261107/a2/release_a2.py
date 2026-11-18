"""Write the A2 reader report and a portable archive with its frozen inputs."""
import csv
from html import escape
import json
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
from check_inputs import HERE, ROOT, A1, digest

def csv_rows(f):
    with (HERE/f).open(encoding='utf-8-sig',newline='') as s:return list(csv.DictReader(s))
def fmt(r,key,n=3): return f'{float(r[key]):.{n}f}'
def table(head,rows):
    h=''.join(f'<th>{escape(k)}</th>' for k in head)
    tr=''.join('<tr>'+''.join(f'<td>{escape(str(v))}</td>' for v in row)+'</tr>' for row in rows)
    return f'<div class="scroll"><table><thead><tr>{h}</tr></thead><tbody>{tr}</tbody></table></div>'
def report():
    z=json.loads((HERE/'a2_verified.json').read_bytes()); ev=csv_rows('event_comparison.csv');tt=csv_rows('trial_summary.csv')
    confirmed=[r for r in ev if r['confirmed']=='1'];early=z['initial'];late=z['later']
    h='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A2 原认证范围内的端点约束删行核查</title><style>body{max-width:1120px;margin:auto;padding:30px min(4vw,50px) 65px;
font:16px/1.72 "Microsoft YaHei",sans-serif;color:#172831}h1{font-size:30px;line-height:1.25}h2{font-size:21px;color:#195169;margin:28px 0 7px}p{margin:8px 0}table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}td,th{text-align:right;padding:5px 8px;border-bottom:1px solid #d2dfe2;white-space:nowrap}th{background:#e6f0f2}td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}a{color:#16576b}small{color:#425965}
@media print{body{padding:0}thead{display:table-header-group}tr{break-inside:avoid}.scroll{overflow:visible}}</style><main>'''
    h+=f'''<h1>A2 原认证范围内的端点约束删行核查</h1>
<p><strong>在原 A1 选定的 {z['trials']} 条 Matched-QP 轨迹中，删除端点行改变了 {z['confirmed']} 个保存点的最优参考指令，其他 {z['exactly_equal']:,} 个点的指令数值完全一致。</strong>
{late['points']} 个变化点分布在飞行过程的 {len(late['trial_ids'])} 条轨迹上。原 QP 当时仍强制备份行和耦合推力限制。所有变化点都能由上端点导数条件解释。</p>
<p>原 500 对试验中，431 条 Matched-QP 轨迹完成；A1 回查的是有上端点正对偶支持的全部 17 条完成轨迹。本次全网格分母 {z['total']:,} 是这 17 条的保存点数。时间步为原输出网格 0.05 s。采样点的参与次数不等于连续作用时长。</p>'''
    def grouprow(name,g):
        return [name,len(g['trial_ids']),g['points'],g['backup_support'],g['J_included'],
                f"{g['d_range'][0]:.4f} 至 {g['d_range'][1]:.4f}",
                f"{g['upper_phi_range'][0]:.3f} 至 {g['upper_phi_range'][1]:.3f}"]
    h+=table(['时段','涉及试验数','确认点','原备份同时支持','A1 联合盒内','指令差 d','删行解代回上端点 φ (N/s)'],
             [grouprow('初始化 0 s',early),grouprow(f"飞行中 {late['time_range'][0]:.2f}–{late['time_range'][1]:.2f} s",late)])
    low=min(float(r['phi16_lower_Nps']) for r in confirmed)
    h+=f'''<p>指令差 d = ‖S⁻¹(u16−u18)‖∞ 为无量纲量。对于全部 {len(confirmed)} 个变化点，16 行解代回下端点导数行的最小余量仍为 +{low:.3f} N/s。由于 16 行解已满足下端点条件，在同一严格凸 QP 中，将下端点行单独加回也保留该最优解。因此可以将这组点的指令差归于上端点行。此处的 φ 是导数条件余量；参考/反馈根的当前端点力间隙以 N 计。</p>'''
    h+='<h2>两个没有原备份正对偶支持的飞行中样本点</h2>'
    highlights=[r for r in confirmed if float(r['time_s'])>0 and r['backup18']=='0']
    h+=table(['试验','t (s)','d','ΔuV (m/s³)','Δuγ (°/s²)','Δuχ (°/s²)','上端点 φ16 (N/s)','删行后备份支持'],
             [[r['id'],fmt(r,'time_s',2),fmt(r,'d_scaled_inf',4),fmt(r,'delta_V_mps3'),fmt(r,'delta_gamma_degps2'),
               fmt(r,'delta_chi_degps2'),fmt(r,'phi16_upper_Nps'),r['backup16']] for r in highlights])
    h+='<p>差值取删行解减原解，三个分量属于参考模型的输入。备份行在两种配置中都存在。表中最后一列表示删行后备份行是否具有正对偶支持。</p>'
    r=next(r for r in confirmed if r['id']=='86' and abs(float(r['time_s'])-68.85)<.0001)
    h+=f'''<p>样本 86 在 68.85 s 也有确认差异，d={fmt(r,'d_scaled_inf',4)}，删行解代回上端点条件为 {fmt(r,'phi16_upper_Nps')} N/s。A1 同时刻直接重建的反馈上端点间隙为 +{fmt(r,'fb_upper_N')} N，未裁剪推力请求为 {fmt(r,'Traw_N',4)} N。在原状态固定的这一拍，反馈使用当前参考状态和实际状态。后续是否改善反馈力裕度、推力裁剪或位置 RMS，需要重新积分两条闭环。</p>'''
    h+='<h2>逐点结果</h2><p>φ18 与 φ16 均针对原上端点行。备份列分别为原解和删行解的备份正对偶支持；J 列沿用 A1 对原轨迹、原保存时刻的条件性判定。</p>'
    for label,sub in [('飞行中',[r for r in confirmed if float(r['time_s'])>0]),('初始化',[r for r in confirmed if float(r['time_s'])==0])]:
        h+=f'<p><strong>{label}</strong></p>'
        h+=table(['试验','t (s)','d','ΔuV (m/s³)','Δuγ (°/s²)','Δuχ (°/s²)','φ18 (N/s)','φ16 (N/s)','备份 18 / 16','J'],
            [[a['id'],fmt(a,'time_s',2),fmt(a,'d_scaled_inf',4),fmt(a,'delta_V_mps3'),fmt(a,'delta_gamma_degps2'),
              fmt(a,'delta_chi_degps2'),f"{float(a['phi18_upper_Nps']):.2g}",fmt(a,'phi16_upper_Nps'),
              a['backup18']+' / '+a['backup16'],a['J_included']] for a in sorted(sub,key=lambda a:(float(a['time_s']),int(a['id'])))])
    h+=f'''<h2>复算方法与验收</h2>
<p>在每一原 MC 状态、时间和外部条件下，W2E 快照重建原请求 u0、原 18 行解和全部余量，并核对 A1 时间轴与 18 行系数。全网格最大请求缩放差、原解缩放差、相对余量重放差分别为 {z['max_original_replay']['max_u0_replay']:.2g}、{z['max_original_replay']['max_u18_replay']:.2g}、{z['max_original_replay']['max_phi_replay']:.2g}。准入及嵌套可行性检查通过 {z['eligible']:,}/{z['total']:,} 点。</p>
<p>仅删 14、15 行，保留原索引 1–13 和 16–18（含推力下界、上界及备份）。两个 QP 共用原目标、尺度 diag(2,8π/180,12π/180) 和可行起点。确认需满足指令差 d&gt;10⁻⁸，并且同一条被删行的归一化违反量&gt;10⁻⁸、原量纲负余量超过该行求解器容差。两解的 KKT 残差≤5×10⁻⁸，保留行归一化违反量≤10⁻⁸，并检查 J16≤J18+10⁻⁸(1+J18)。这里的 J16、J18 是优化目标，与 A1 的反馈充分盒 J 分开记录。</p>
<p>对原端点支持点、确认候选及数值边缘或判据分歧点触发枚举；数值边缘采用相应阈值的 0.5–2 倍。独立复核包含 MATLAB 全活跃集 QR 枚举和 Python SVD 枚举（支持规模 0–3），两种配置均求解。共枚举 {z['enumerated']} 个点。两项枚举与生产解的最大缩放指令差为 {z['max_solver_comparison']:.4g}。Python 另在全网格以对偶、原始行和状态矩阵重算 KKT，18/16 行的最差归一化残差为 {z['max_independent_kkt'][0]:.4g}/{z['max_independent_kkt'][1]:.4g}。排除或判据分歧 {z['exceptions']} 点。MAT 文件保存两解和枚举解，CSV 保留有符号余量和门槛。</p>'''
    ids='、'.join(map(str,z['a3_ids'])); first={int(a['id']):a['first_confirmed_s'] for a in tt}
    h+=f'''<h2>论文与后续实验</h2><p>现有证据支持“上端点导数条件在原备份和推力限制仍生效时，可独立改变 Matched-QP 的瞬时选择”。这些保存点处在原备份的数值准入门槛内。A1 在所有 21 点均观察到原反馈反演有根。A2 固定状态，无法给出删行闭环的后续反馈裕度、实际执行器状态或跟踪代价，也未检验 Scalar 或相对等预留攻角保护的性能优势。</p>
<p>A3 的进入条件已达到。按方案在全部确认编号中升序取前三例，编号为 {ids}，原轨迹首个确认时刻依次为 {", ".join(f"{float(first[i]):.2f}" for i in z['a3_ids'])} s。A2 交付只包含代数回放。</p>
<p><a href="trial_summary.csv">17 例汇总</a>　<a href="event_comparison.csv">事件明细及 A1 对照</a>　<a href="all_samples.csv">全部保存点</a>　<a href="enumeration_checks.csv">独立枚举</a>　<a href="exceptions.csv">异常账本</a>　<a href="README_复算说明.txt">数据与复算说明</a></p></main></html>'''
    (HERE/'A2_同状态删行核查.html').write_text(h,encoding='utf-8')

def release():
    report()
    text="""Aerospace 大修 A2 同状态删行核查

打开 A2_同状态删行核查.html。trial_summary.csv 为 17 个冻结样本的逐例验收，event_comparison.csv 为事件级数据，all_samples.csv 为全网格。MATLAB 索引从 1 开始。枚举表含 QR 及独立 SVD 核查。exceptions.csv 若只有表头，表示该批次无排除/分歧点。A3 编号按原方案确定，仅用于交接。

结果字段
u0/u18/u16 是原请求、18 行最优指令和删行后的 16 行最优指令（三列依次为速度、爬升角、航向的二阶导数），MAT 中单位 m/s³、rad/s²、rad/s²。delta 均取 u16-u18，CSV 为阅读方便将两个角分量转为 °/s²。d_scaled_inf 为缩放差的无穷范数；phi18/phi16 是把相应指令代回原端点导数行的值（N/s）。事件级 ref_upper_N、fb_upper_N 及 Traw_N 来自 A1，只描述该拍原状态。
runs/id_XXX.mat 保存全时序、数值结果矩阵 M 及列名 fields、状态字符串 status、三类指令、删行对偶及触发点的两种 MATLAB 枚举解。未触发枚举的数值列用 NaN。objective18/16 为对应同状态优化目标。

复算
解压后保持目录结构。原研究参数、A1 轨迹与系数、W2E 和原机模源码均随包提供。使用 Python 3.11+、NumPy、SciPy，运行 a2 目录下的 python check_inputs.py；在 MATLAB R2018a 的 a2 目录运行 run_a2，再执行 python verify_a2.py。前者重新构造并求解两套 QP，后者独立重算数值判据、穷举候选及阈值分歧点并输出 CSV/JSON。如仅需验证已有数值结果，直接运行 python verify_a2.py。release_a2.py 生成阅读报告及本地归档。
全部计算只重建冻结状态的代数接口。重跑时在独立副本内操作，以便保留初次结果。input_hashes.json 锁定本次使用的 A1 输入与原代码，release_manifest.json 提供本包文件大小与 SHA-256。
"""
    (HERE/'README_复算说明.txt').write_text(text,encoding='utf-8-sig')
    items=set(HERE.glob('*.py'))|set(HERE.glob('*.m'))|set(HERE.glob('*.csv'))|set(HERE.glob('*.txt'))|set(HERE.glob('*.html'))|set(HERE.glob('*.log'))
    items.update([HERE/'input_hashes.json',HERE/'a2_verified.json',A1/'release_manifest.json'])
    items.update((HERE/'runs').glob('*.mat'));items.update((HERE/'identity').glob('*.csv'))
    items.update(ROOT/k for k in json.loads((HERE/'input_hashes.json').read_bytes()))
    items.add(ROOT/'revision_stage1_20261107/baseline/code/P1_Aerospace_W2B_R2018a_v1_1/w2b_admission.m')
    releasefile=HERE/'release_manifest.json'
    manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(items)}
    releasefile.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    items.add(releasefile)
    final=ROOT/'Aerospace_A2同状态删行核查材料包_v1_0.zip';tmp=final.with_suffix('.partial.zip')
    with ZipFile(tmp,'w',compression=ZIP_DEFLATED,compresslevel=3) as z:
        for f in sorted(items):z.write(f,f.relative_to(ROOT).as_posix())
    with ZipFile(tmp) as z:
        assert set(z.namelist())=={p.relative_to(ROOT).as_posix() for p in items}
        for k,v in manifest.items():
            with z.open(k) as f:assert hashlib.file_digest(f,'sha256').hexdigest()==v['sha256']
    tmp.replace(final);print(len(items),'files;',final.stat().st_size,'bytes;',digest(final))
if __name__=='__main__':release()
