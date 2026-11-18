"""Generate the aligned tables, Chinese research note and integrity manifest."""
import html,json,zipfile
from pathlib import Path
import pandas as pd
import numpy as np
from freeze import HERE,RESULT,ROOT,lock,digest

def table(headers,rows):
    def cells(tag,vals):return '<tr>'+''.join(f'<{tag}>{html.escape(str(v))}</{tag}>' for v in vals)+'</tr>'
    return '<table><thead>'+cells('th',headers)+'</thead><tbody>'+''.join(cells('td',r) for r in rows)+'</tbody></table>'

def run():
    lock()
    def read(n):return pd.read_csv(RESULT/n,float_precision='round_trip')
    a=read('actual_feedback.csv');lp=read('historical_lp.csv');f=read('feedback_boxes.csv');sec=read('feedback_section_probe.csv')
    rs=read('recovery_rows.csv');corner=rs[rs.case.str.startswith('section')]
    r=rs[~rs.case.str.startswith('section')].copy()
    r['key']=r.time_s.round(10)
    # Use physical seconds, not index. Never interpolate one arm to another terminal event.
    for g,prefix,drop in [(lp,'LP',['VB_over_rho']), (a,'observed',[])]:
        g=g.drop(columns=drop).copy();g['key']=g.time_s.round(10);g=g.drop(columns='time_s')
        g=g.rename(columns={k:prefix+'_'+k for k in g if k not in ['case','key']})
        r=r.merge(g,on=['case','key'],how='left',validate='one_to_one')
    for b,g in f.groupby('box'):
        g=g.copy();g['key']=g.time_s.round(10);g=g.drop(columns=['time_s','box'])
        g=g.rename(columns={k:b+'_'+k for k in g if k not in ['case','key']})
        r=r.merge(g,on=['case','key'],how='left',validate='one_to_one')
    r=r.drop(columns='key').sort_values(['time_s','case'])
    r.to_csv(RESULT/'time_aligned.csv',index=False,encoding='utf-8-sig')
    points=[5.55,5.6,7.75,7.8,12,14,14.60,14.65,14.70,15.55,16,16.8,18.7,20,24.8,35,55]
    keep=r.time_s.round(8).isin(points)|(r.time_s.between(14.650001,14.653))
    cols=['time_s','case','Tref_N','hNplus','LP_backup_negative','LP_lp_margin_scaled','LP_uV','LP_uGamma','LP_uChi','observed_ref_upper_N','observed_fb_upper_N','observed_loss_lower_N','observed_loss_upper_N','directional_lower_loss_sup_ub_N','directional_upper_loss_sup_ub_N','directional_upper_fb_inf_lb_N','directional_pay_10N','directional_E','uniform_upper_loss_sup_ub_N','uniform_E']
    r.loc[keep,cols].to_csv(RESULT/'key_times.csv',index=False,encoding='utf-8-sig')
    combined=corner.merge(sec,on=['case','sample_index'],suffixes=('','_iv'),validate='one_to_many')
    assert np.max(abs(combined.time_s-combined.time_s_iv))<1e-12
    combined.drop(columns='time_s_iv').to_csv(RESULT/'sections_aligned.csv',index=False,encoding='utf-8-sig')
    checker=json.loads((RESULT/'verification.json').read_bytes());assert checker['status']=='passed'
    rec=json.loads((RESULT/'recovery_summary.json').read_bytes())
    bad=json.loads((RESULT/'candidate_witnesses.json').read_bytes())
    candidate=[];fb_rows=[];counts=[];summary={}
    for d in [18,22,26]:
        for rho in [30,10,3]:
            k=f'section{d}_{rho:03}'; g=sec[(sec.case==k)&(sec.box=='directional')];v=rec[k]
            s=[d,f'{rho/1000:.3f}',f"{v['positive_all']}/{v['n']}",f'{(g.pay_10N=="yes").sum()}/{len(g)}',f'{g.upper_loss_sup_ub_N.max():.4f}']
            if k in bad:
                q=bad[k];extra='，'.join(('推力下限' if c=='T_min' else '上端点预留')+f" {z['reserve_margin_N'][1]:.4f} N（{z['time_s']:g} s）" for c,z in q.items())
                decision='已排除。'+extra
            elif (g.pay_10N=='counterexample').any():decision='恢复采样通过；10 N 有反例，所测点仍括根。'
            else:decision='首选进入区域验证' if k=='section22_010' else '保留，依原顺序候补'
            candidate.append(s+[decision]);summary[k]=dict(recovery=v['positive_all']==v['n'],budget_counterexamples=int((g.pay_10N=='counterexample').sum()),upper_loss_bound_N=g.upper_loss_sup_ub_N.max(),samples=len(g),recommendation=decision)
    for nm in ['N10','T','center18','center22','center26']:
        for box in ['directional','uniform']:
            g=f[(f.case==nm)&(f.box==box)&(f.time_s.between(14,20))]
            fb_rows.append([nm,{'directional':'分轴','uniform':'统一'}[box],len(g),f'{g.lower_loss_sup_ub_N.max():.4f}',f'{g.upper_loss_sup_ub_N.max():.4f}',f'{g.upper_fb_inf_lb_N.min():.4f}',(g.pay_10N=='yes').sum(),(g.E=='yes').sum()])
    for case,g in lp.groupby('case',sort=False):
        for begin,end in [(0,35),(3,16),(14,20)]:
            sub=g[g.time_s.between(begin,end)];counts.append([case,f'{begin}–{end}',len(sub),int((sub.VB_over_rho<=1).sum()),int(sub.backup_negative.isna().sum()),int(sub.lp_feasible.sum()),f'{sub.lp_margin_scaled.min():.6f}'])
    outcome={'recommended':{'decel_center_s':22,'rho':.01,'center_target_mps':31.25},'candidates':summary,'status':'Section 2 screens completed; Section 3 regional certificate pending','confirmation_runs':0}
    (RESULT/'decision.json').write_text(json.dumps(outcome,indent=2,ensure_ascii=False),encoding='utf8')
    css='body{font:16px/1.65 "Source Sans 3","Microsoft YaHei",sans-serif;color:#152932;max-width:1080px;margin:32px auto;padding:0 26px}h1{font-size:31px}h2{font-size:22px;margin:32px 0 8px}p{margin:9px 0}table{border-collapse:collapse;width:100%;margin:12px 0;font-variant-numeric:tabular-nums;font-size:14px}th,td{border-bottom:1px solid #ccd5d8;padding:6px 7px;text-align:left}th{background:#e7eff1}a{color:#075482}.decision{border-left:4px solid #28665e;padding:3px 14px;background:#eef6f3}.meta{color:#46616b}@media print{body{font-size:10.5pt}h2{break-after:avoid}tr{break-inside:avoid}table{font-size:9pt}}'
    body=['<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><title>有限时长减速机动 · 两项离线筛查</title><style>'+css+'</style><main>',
      '<h1>有限时长减速机动的两项离线筛查</h1><p class="meta">2026 年 11 月 8 日　｜　Aerospace 修订研究　｜　协议第 2 节</p>',
      '<p class="decision"><b>进入区域验证的首个候选为 22 s 减速中心、ρ=0.01。</b>该中心从 38.1 减至 31.25 m/s，35–55 s 的恢复尾段参加筛查。22 s 是备份中心的时长。原拟直线地速请求仍用 8 或 10 s 减速。此次计算全部基于历史记录和解析参考，新增闭环运行数为零。</p>',
      '<h2>原备份何时失去作用</h2><p>H31 的 N10 轨迹从初值到 5.55 s 仍在原椭球。5.60 s 是首个域外保存点，直至 7.80 s 才首次出现原备份输入的负导数行。该点的推力上限行残差为 −0.579048 N/s，实际参考推力仍在 3–57 N 内。同点换用线性规划给出的输入，推力上限行残差为 +9.910327 N/s，全部 17 行同时满足。16.80 s 原备份的上端点行残差为 −81.983626 N/s，此时历史轨迹已经开始转向。</p>',
      table(['历史记录','统计窗 / s','保存点','原椭球内','原备份 17 行通过','另求可行输入','最小归一余量'],counts),
      '<p>线性规划最大化 17 行的共同归一余量，归一分母为 max(1,‖Gⱼ diag(scale)‖₂)，三项输入按原 QP 尺度缩放。未另设加速度限幅。每个可行解均代回 MATLAB 的原始 c+Gu，原单位残差可逐行复算。此表检验当前点存在输入；单点解未规定下一时刻的状态，也不能修复已违反的状态约束。T 在 14.652183 s 终止，统计窗按已有记录截断。</p>',
      '<h2>哪些候选在正宽度截面上仍成立</h2><p>对每个中心沿 0–55 s 每 0.05 s 计算一次（1,101 点），所有中心的状态余量和运行导数余量均为正。在预定相位按 Cholesky 坐标的五条轴取正负两侧，偏移达半径的 97%。全程合格列要求同点满足所有状态、初始化、根正则及 17 条备份导数条件。关键窗两列只统计 14–20 s 已取的非零偏移点，并在各点检验整个分轴误差盒。</p>',
      table(['中心 / s','ρ','全程截面合格点','10 N 可支付点','上端损失上界 / N','判定或严格负见证'],candidate),
      '<p>例如，18 s、ρ=0.003 在 12 s 的一个管内点只需 2.580017 N 参考推力，低于 3 N 设计预留。22 s、ρ=0.030 和 26 s、ρ=0.030 的低速截面都出现上端点预留短缺 3.752929 N。各负见证重新按解析中心和原二进制参数外向计算，包含域成员、攻角根括区、斜率和推力区间。首选候选所测截面最低参考推力为 3.869234 N，最小上端点状态余量为 3.402151 N（已扣除 10 N）。</p>',
      '<h2>反馈端的力预算</h2><p>H31/N10 的 35 s 全迹三轴最大误差为 0.122151 m/s、0.419992°、0.329437°。按协议向上圆整得分轴盒 ±(0.125 m/s, 0.45°, 0.35°)。统一盒取相同角度尺度 38.1 m/s，容纳上述三轴的最小速度半宽约为 0.29924 m/s。每个时间点计算有符号损失 ℓ=g参考−g反馈。两端损失上界均小于 10 N 记为可支付。反馈根栏要求两个原始反馈力差及力映射斜率的下界均为正。</p>',
      table(['参考点组','误差盒','关键窗点数','下端损失上界 / N','上端损失上界 / N','上端反馈下界 / N','可支付','可括根'],fb_rows),
      '<p>在 18.70 s 的原转向记录，N10 上端参考裕度为 10.144552 N，实际消耗 6.707683 N，实际剩余 3.436869 N。对此参考点的整个分轴盒，最坏上端损失约为 9.229211 N，反馈端仍留约 0.915340 N。整个关键窗的 121 个 N10 保存点中，20 s 对分轴盒最苛刻，反馈上端下界约为 0.249870 N。统一盒在同窗 121 个点均可找到损失超过 10 N 的允许误差，其中 60 个点存在失去反馈括根的误差。这里描述的是盒内允许误差，不能解读为沿实测轨迹发生了 60 次退出。</p>',
      '<p>H31/T 在 14.652183 s 触及反馈上端边界：参考端剩 5.844171 N，实际损失 5.844171 N。此时参考推力为 7.609753 N，原始反馈推力为 7.607063 N。该点说明推力合规仍能耗尽反馈反演裕度。表中 T 的“可支付”只表示给足 10 N 时误差损失落在预算内；T 没有施加这项预留。两条历史对照的期望位置与速度在转向前使用同一地面请求；已核对 H31 和直线请求在 0–16 s 的 321 个网格输入逐值相同。直线闭环的 16 s 后状态尚未生成。</p>',
      '<h2>放行范围</h2><p>22 s、ρ=0.01 的关键窗内，方向盒在 121 个中心点及 70 个预定截面点均通过。截面点上端损失上界约为 9.795603 N，距 10 N 约 0.2044 N；同批点使用统一盒存在损失约 12.3875 N 的反例，但额外名义力裕度足以保持括根。在整个连续时间与完整五维截面上，结论仍待第 3 节的区域计算。优先检查上端力损失、3 N 推力预留和上端点预留，再完成全部运行行、含显式时间导数的第 18 行及终端包含性。G/J 联合条件、实际误差持续驻盒、C1–C3 的工程效益均留待相应阶段。</p>',
      '<h2>算法与复核</h2><p>原备份使用冻结配置 ρ=0.10。反馈误差盒的区间计算采用 40 位十进精度；对完整攻角区间分 8 段检验正斜率。先按航向误差对合力的单调性找出最不利面，再在速度与爬升角偏导严格同号时缩到相应端面，余下维度保留区间。可支付判断使用整个盒的向外上界，严格反例使用允许顶点的损失下界。独立浮点公式重查了 105,786 个参考状态与误差组合，最大超出区间量约 3×10⁻¹³ N（数值复核容差 10⁻⁹ N）。连续时间、完整移动管与实际六自由度安全均未在本阶段认证。</p>',
      '<p><a href="results/key_times.csv">关键时刻对齐表</a>　<a href="results/time_aligned.csv">全部轨迹与中心对齐表</a>　<a href="results/sections_aligned.csv">关键窗截面对齐表</a>　<a href="results/candidate_witnesses.json">严格负见证</a>　<a href="results/row_minima.csv">逐行最差位置</a>　<a href="README_复算.txt">复算说明</a></p></main></html>']
    historic=r[r.case.isin(['N10','T'])]
    events=historic[((historic.case=='N10')&historic.time_s.round(8).isin([5.60,7.80,14.65,14.70,16.80,18.70,20]))|
                    ((historic.case=='T')&(historic.time_s.between(14.649999,14.653)))]
    key_rows=[]
    for _,v in events.iterrows():
        negatives=v.LP_backup_negative
        key_rows.append([f'{v.time_s:.6f}'.rstrip('0').rstrip('.'),v['case'],
                         negatives if pd.notna(negatives) else '均满足',f'{v.LP_lp_margin_scaled:.4f}',
                         f'{v.observed_loss_lower_N:.4f}',f'{v.observed_loss_upper_N:.4f}',
                         f'{v.observed_fb_upper_N:.4f}',str(v.directional_pay_10N).replace('yes','可支付')])
    insert_at=next(i for i,t in enumerate(body) if t.startswith('<h2>哪些候选'))
    body[insert_at:insert_at]=['<h2>历史关键时刻对齐表</h2><p>导数负行针对假设施加原备份的瞬间。两端力损失及上端剩余量来自实际执行的历史控制。最后一列单独检查整个分轴盒的误差损失。事件时刻只列原始记录已有的方案。</p>',
        table(['时间 / s','方案','原备份负行','另求输入最小归一余量','实测下端损失 / N','实测上端损失 / N','实测上端剩余 / N','10 N 预算'],key_rows)]
    (HERE/'两项离线筛查报告.html').write_text('\n'.join(body),encoding='utf8')
    parts=['freeze.py','protocol_lock.json','export_affine.m','screen_recovery.py','screen_feedback.py','check_sections.py','verify.py','assemble.py','README_复算.txt','两项离线筛查报告.html']
    files=[HERE/p for p in parts]+sorted(p for p in RESULT.iterdir() if p.is_file())
    manifest={str(p.relative_to(HERE)):digest(p) for p in files}
    fp=HERE/'release_manifest.json';fp.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf8')
    pack=ROOT/'Aerospace_有限时长减速机动_第2节离线筛查_v1_0.zip'
    with zipfile.ZipFile(pack,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in [*files,fp]:z.write(p,str(p.relative_to(HERE)))
    with zipfile.ZipFile(pack) as z:
        assert z.testzip() is None and len(z.namelist())==len(files)+1
    print(pack.name,pack.stat().st_size)
if __name__=='__main__':run()
