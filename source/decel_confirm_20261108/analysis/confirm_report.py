"""Assemble cross-case decisions and an auditable report from the checked traces."""
from html import escape
import json, zipfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from confirm_metrics import H,BASE,OUT,ROOT,main as evaluate,digest

REPORT='C1–C3_直线减速确认报告.html'
def table(head,rows):
    def tr(tag,row):
        return '<tr>'+''.join(f'<{tag}>{escape(str(cell))}</{tag}>' for cell in row)+'</tr>'
    return '<table><thead>'+tr('th',head)+'</thead><tbody>'+''.join(tr('td',row) for row in rows)+'</tbody></table>'
def fmt(v,n=4):
    return f'{v:.{n}f}' if pd.notna(v) else '—'
def read(name):
    return pd.read_csv(OUT/(name+'.csv'),dtype={'step':str})
def figure(tr):
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(2,2,figsize=(10.6,7.3),layout='constrained')
    for ax,sc,txt in zip(axes.flat,['C1','C2','C3'],['31 m/s，10 s，静风','31.5 m/s，10 s，静风','31.5 m/s，10 s，弱时变风']):
        f=tr[sc,'F','005'];o=tr[sc,'O','005']
        ax.axvspan(3,16,color='#edf3f5');ax.plot(f.time,f.vrequest,'--',color='#364456',lw=1.3,label='请求地速')
        ax.plot(o.time,o.vground,color='#ad653a',lw=1.7,label='O 实际地速')
        ax.plot(f.time,f.vground,color='#00776c',lw=1.9,label='F 实际地速')
        ax.plot(o.time.iloc[-1],o.vground.iloc[-1],marker='x',color='#ad653a',ms=6)
        if sc=='C3':ax.plot(f.time,f.airSpeed,':',color='#704fa1',lw=1.6,label='F 实际空速')
        ax.set(title=f'{sc}   {txt}',xlim=(0,35),ylim=(29.8,39),ylabel='速度 (m/s)',xlabel='时间 (s)')
        ax.legend(frameon=False,fontsize=8)
    ax=axes.flat[-1]
    for arm,c,ls in [('T','#ad653a','-'),('N10','#00776c','-'),('Am','#584c99','--')]:
        z=tr['C1',arm,'005'];ax.plot(z.time,z.fbHi,ls,color=c,lw=1.6,label=arm)
        if arm=='T':ax.plot(z.time.iloc[-1],z.fbHi.iloc[-1],'o',color=c,ms=5)
    ax.axhline(0,color='#364456',lw=.7);ax.axvspan(14,20,color='#edf3f5',zorder=-1)
    ax.set(title='C1   去除备份行后的反馈上端力差',xlim=(13,20.1),ylim=(-2,62),ylabel='反馈上端力差 (N)',xlabel='时间 (s)')
    ax.legend(frameon=False,fontsize=8)
    for ax in axes.flat:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15)
    fig.savefig(OUT/'confirmation.png',dpi=160);plt.close(fig)

def assemble():
    metrics,tr=evaluate(); figure(tr)
    cases=['C1','C2','C3']; idx=metrics.set_index(['case','arm','step']).sort_index()
    def row(sc,arm):return idx.loc[sc,arm,'005']
    c=read('contrasts').set_index(['case','name'])
    cover=read('coverage'); ab=read('matched_replay');free=read('free_replay'); common=read('common_windows')
    changes=read('cross_case_F').set_index(['from_case','to_case','metric'])
    def difference(a,b,metric):return changes.loc[a,b,metric]
    evt=read('events'); ac=read('actuators')
    decision={}
    for sc in cases:
        def get_both(arm):return idx.loc[sc,arm]
        f=get_both('F');o=get_both('O')
        dec=dict(speed_scale_pass=bool(c.loc[sc,'O_minus_F_speed'].meets_scale),both_F_complete=bool((f.cause=='completed').all()),
                 out_of_old_low_plateau_both=bool(((f.low_plateau_samples==201)&(f.plateau_samples==201)).all()),O_events=o.cause.unique().tolist())
        for alt in ['D','A']:
            altrow=get_both(alt);same=ab[(ab['case']==sc)&(ab.substitute==alt)]
            changed=bool((same.changed>0).all())
            hit=all(altrow.loc[step].cause in ['feedback_lower_root','feedback_upper_root'] and f.loc[step].cause=='completed' and min(f.loc[step].fbHi_min,f.loc[step].fbLo_min)>0 for step in ['010','005'])
            dec['independent_over_'+alt]=changed and (hit or bool(c.loc[sc,'F_minus_'+alt+'_upper'].meets_scale))
        decision[sc]=dec
    strong=all(v['speed_scale_pass'] and v['both_F_complete'] and v['out_of_old_low_plateau_both'] and v['independent_over_D'] and v['independent_over_A'] for v in decision.values())
    assert not strong,'Strong outcome requires the conditional Scalar runs before release.'
    et=evt.query("case=='C1' and arm=='T'").set_index('step')
    te=et.loc['005'];nt=row('C1','N10');am=row('C1','Am')
    mech=bool((et.cause=='feedback_upper_root').all() and abs(et.loc['010','event_s']-te.event_s)<.05 and
        (idx.loc['C1','N10'].cause=='completed').all() and (idx.loc['C1','N10'].fbHi_min>0).all() and (free.changed_T>0).all())
    full=dict(unique_case_arms=len(metrics)//2,integrations=len(metrics),completed35=int((metrics.cause=='completed').sum()),results=decision,
        C1_no_tube_N10_vs_T_reproduced=mech,C1_no_tube_N10_vs_Am_scale=bool(c.loc['C1','N10_minus_Am_upper'].meets_scale),
        C1_C2_governed_plateau_delta=float(difference('C1','C2','vref_ground_25_35').fine),
        C1_C2_plateau_resolved=bool(difference('C1','C2','vref_ground_25_35').resolved),
        strong_gate=strong,optional_scalar_integrations=0)
    (OUT/'decision.json').write_text(json.dumps(full,ensure_ascii=False,indent=2),encoding='utf8')
    # The prose below describes the measured branch of the prespecified comparisons.
    assert mech and all(v['speed_scale_pass'] and v['out_of_old_low_plateau_both'] for v in decision.values())
    assert (ab.changed==0).all() and (ab.max_abs_36state_diff==0).all() and ab.same_grid.all()
    assert not full['C1_no_tube_N10_vs_Am_scale'] and not full['C1_C2_plateau_resolved']
    critical=cover.query("window=='critical' and arm=='F'")
    assert ((critical.n==121)&(critical.EG==121)&(critical.J==121)).all()
    sums=metrics[metrics.arm=='F']
    assert (sums.interventions==sums.backup_support_on_interventions).all() and (sums.endpoint_support==0).all()
    df=pd.read_csv(BASE/'results/run_metrics.csv',dtype={'step':str})
    hist=df[df.step=='005'].set_index('arm')
    base_delta=hist.loc['O'].MAE_3_16-hist.loc['F'].MAE_3_16
    def p(content):return '<p>'+content+'</p>'
    doc=['''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>C1–C3 减速机动确认</title>
<style>body{max-width:1050px;margin:3em auto;padding:0 1.3em;color:#18313b;background:#fff;font:16px/1.68 "Microsoft YaHei",sans-serif}h1,h2{line-height:1.35}h1{font-size:2em}h2{margin-top:1.6em;font-size:1.25em}table{border-collapse:collapse;width:100%;font-size:.9em;margin:.8em 0;display:block;overflow-x:auto}td,th{padding:.4em .65em;border-bottom:1px solid #d7e0e2;text-align:right;white-space:nowrap}th{background:#edf3f5}td:first-child,th:first-child{text-align:left}img{max-width:100%}a{color:#00695c}.lead{border-left:4px solid #00776c;padding-left:1em}</style>
<h1>有限时长直线减速机动 · C1–C3 确认结果</h1>''',
    '<p class="lead">三个确认工况均复现了移动备份相对原静态备份的减速改善。相同移动管中，完整端点、删端点和等预留攻角三组产生相同的指令及轨迹。取消备份行后，C1 重现了仅限推力时的反馈上端括根事件，完整端点与等预留攻角均避开该事件。</p>',
    '<h2>工况与计量</h2>',p('三个确认工况均从 38.1 m/s 配平开始，3 s 起按五次光滑曲线在 10 s 内减速，保持至 35 s。C1、C2 依次请求 31、31.5 m/s，采用静风。C3 沿用 C2 请求，加入预定的弱时变风及其解析导数。全程关闭噪声。原静态备份使用 ρ=0.10；三组移动备份共享 22 s 时间表、31.25 m/s 平台中心及 ρ=0.01。F/O 使用相同完整端点及同一 QP，D 删除两条端点行并保留推力限制，A 采用等 10 N 的移动攻角界。'),
    p(f'共 15 个不同的工况–配置组合，各以最大积分步长 0.01、0.005 s 独立运行，共 30 条轨迹，{full["completed35"]} 条到达 35 s。下表引用细步长；汇总文件同时保留粗步长。0.05 s 采样，以梯形法计算时间平均。原静态组在原始指令坐标门槛处提前退出；位置信息另外在配对共同网格上比较。'),
    '<h2>改变减速时长、目标和风后的地速表现</h2>']
    mainrows=[['D0 · 31 / 8 / 静风',fmt(hist.loc['O'].MAE_3_16),fmt(hist.loc['F'].MAE_3_16),fmt(base_delta),fmt(hist.loc['O'].end_time,3),fmt(hist.loc['F'].vground_25_35),fmt(hist.loc['F'].governed_25_35)]]
    for sc,target in zip(cases,['31 / 10 / 静风','31.5 / 10 / 静风','31.5 / 10 / 弱风']):
        f=row(sc,'F');o=row(sc,'O')
        mainrows.append([sc+' · '+target,fmt(o.MAE_3_16),fmt(f.MAE_3_16),fmt(c.loc[sc,'O_minus_F_speed'].fine),fmt(o.end_time,3),fmt(f.vground_25_35),fmt(f.vref_ground_25_35)])
    doc.append(table(['目标地速 / 时长 / 风','O 地速 MAE','F 地速 MAE','改善 O−F','O 退出 (s)','F 平台实际地速','F 平台治理地速'],mainrows))
    drift=c.xs('O_minus_F_speed',level='name').drift_gate.max()
    doc.append(p(f'MAE 窗为 3–16 s，平台窗为 25–35 s，速度单位均为 m/s。三个改善量均超过预定研究尺度 0.606 m/s 和数值漂移门槛 {drift:.5f} m/s。O 的三次退出均为原始速度指令模触及 5 m/s；参考与反馈的上下端力差在该边界上均为正。F 完成三个工况，三组的低速平台各有 201/201 个保存点位于原认证域的低速侧。'))
    r=difference('C1','C2','vref_ground_25_35'); g=difference('C1','C2','vground_25_35');f3=row('C3','F')
    doc.append(p(f'目标从 31 提高到 31.5 m/s 后，F 的平台治理地速均值只改变 {r.fine:+.8f} m/s，实际地速均值仅改变 {g.fine:+.8f} m/s，均低于 {r.drift_gate:.5f} m/s 数值门槛。静风下名义参考平台的空速投影约为 [31.05837, 31.44163] m/s，两个请求落在该投影的两侧。低速段体现了这张固定时间表允许实现的减速；该对照缺乏对两个相邻请求终值的分辨能力。'))
    doc.append(p(f'从 C2 加入弱风得到 C3，F 的平台实际地速改变 {difference("C2","C3","vground_25_35").fine:+.4f} m/s。C3 实际地速与空速均值分别为 {f3.vground_25_35:.4f} 与 {f3.airSpeed_25_35:.4f} m/s，治理参考用估计风换算的地速为 {f3.vref_ground_25_35:.4f} m/s，参考空速为 {f3.refAirSpeed_25_35:.4f} m/s。C3 的真风模在保存网格上最大 {f3.wind_norm_max:.4f} m/s。'))
    doc.extend(['<img src="results/confirmation.png" alt="三个确认工况中原静态和新移动备份的实际地速及 C1 无管组反馈上端力差">',p('地速图中叉号为原静态备份的最后保存点。力差图中圆点为 T 的上端括根边界，虚线为等预留攻角组。'),'<h2>端点行是否独立改变闭环</h2>'])
    aa=[]
    for sc in cases:
        f=row(sc,'F');x=ab[(ab['case']==sc)&(ab.step=='005')].set_index('substitute')
        aa.append([sc,fmt(f.min_fb_hi_window),fmt(f.min_fb_lo_window),f'{int(x.loc["D"].changed)}/{int(x.loc["D"].n)}',f'{int(x.loc["A"].changed)}/{int(x.loc["A"].n)}',f'{int(f.backup_support_on_interventions)}/{int(f.interventions)}',fmt(f.fbHi_min),fmt(f.fbHi_min_time,2)])
    doc.append(table(['工况','14–20 s 上端最小 (N)','下端最小 (N)','改用 D 的指令变化','改用 A 的指令变化','备份支持 / 干预点','全程上端最小 (N)','发生时间 (s)'],aa))
    doc.append(p('两步长、三个工况共六组同状态重解，按归一化输入无穷范数大于 10⁻⁷ 计数。F/D/A 分别积分的 701×36 个状态在各自相同步长下也逐值相同。两端点行都没有获得正对偶支持。同管比较中，两端反馈力差的差值均为 0 N。备份输入、QP 残差和每行的重解结果均可追溯至原始记录。'))
    z=free[free.step=='005'].iloc[0];k=common.query("case=='C1' and step=='005' and arm_a=='N10' and arm_b=='T'").iloc[0]
    doc.append(p(f'C1 独立无管对照中，T 在 {te.event_s:.6f} s 触及反馈上端。此时参考前馈推力、原始反馈推力请求及执行推力依次为 {te.Tref:.4f}、{te.Traw:.4f}、{te.Tact:.4f} N，推力界内出现了反馈可逆区间耗尽。在最后共同网格点 {k.critical_to:.2f} s，T 与 N10 的上端力差分别为 {k.last_upper_b:.4f} 与 {k.last_upper_a:.4f} N。N10 状态上改用 T 时，{int(z.changed_T)}/{int(z.n)} 个点的输入发生可辨变化，始于 {z.first_T:.2f} s，违反被删上端导数行的点有 {int(z.violated_hi_by_T)} 个。'))
    doc.append(p(f'N10、Am 均完成 C1。关键窗最小上端力差分别为 {nt.min_fb_hi_window:.6f}、{am.min_fb_hi_window:.6f} N，差 {nt.min_fb_hi_window-am.min_fb_hi_window:.6f} N；D0 的相应差为 0.001973 N。C1 的差通过双步长漂移门槛，但低于预定 1 N。C2 和 C3 按计划只运行四个主配置，独立无管归因限于 D0、C1。'))
    doc.append('<h2>反馈证书覆盖、位置与执行器</h2>')
    cc=critical[critical.step=='005'].set_index('case');outside=cover.query("case=='C3' and arm=='F' and step=='005' and window=='saved'").iloc[0]
    doc.append(table(['工况','E/G 采样','J 采样','最大 |eV| (m/s)','最大 |eγ| (°)','最大 |eχ| (°)','最大姿态角和 (°)','最大上端损失 (N)'],[[sc,f'{int(cc.loc[sc].EG)}/121',f'{int(cc.loc[sc].J)}/121',*[fmt(cc.loc[sc][key],5) for key in ['max_abs_eV','max_abs_eg_deg','max_abs_ec_deg','theta_sum_max_deg','upper_loss_max_N']]] for sc in cases]))
    doc.append(p(f'沿用第 3 节的分轴盒 |eV|≤0.125 m/s、|eγ|≤0.45°、|eχ|≤0.35°，J 再要求 θF+θR≤0.70°。三组 F 在 14–20 s 的两个步长下均为 121/121，从而在实际保存点上可以调用区域证书。C3 全程仅 {int(outside.airBox)}/{int(outside.n)} 个保存点在分轴盒内，首次出盒 {outside.first_out_box:.2f} s。全程最大空速误差为 {outside.max_abs_eV:.4f} m/s。证书给出满足前提时的连续模型结论，轨迹的入盒核查分辨率为 0.05 s。'))
    comp=common.query("step=='005' and arm_a=='F' and arm_b=='O'").set_index('case')
    doc.append(table(['工况','完整 F 位置 RMS (m)','F 峰值 (m)','共同末点 (s)','共同窗 F / O 位置 RMS (m)','F 实际攻角 (°)','原始反馈推力 (N)','执行推力 (N)'],[[sc,fmt((f:=row(sc,'F')).position_RMS35,3),fmt(f.peak_pos_until_stop,3),fmt(comp.loc[sc].common_to,2),f'{comp.loc[sc].pos_rms_a:.3f} / {comp.loc[sc].pos_rms_b:.3f}',f'{f.alphaActual_deg_min:.3f}–{f.alphaActual_deg_max:.3f}',f'{f.Traw_min:.3f}–{f.Traw_max:.3f}',f'{f.Tact_min:.3f}–{f.Tact_max:.3f}'] for sc in cases]))
    doc.append(p(f'C1 有相同的完整端点保护可供比较。移动管 F 的全程位置 RMS 为 {row("C1","F").position_RMS35:.3f} m，无管 N10 为 {nt.position_RMS35:.3f} m，对应 3–16 s 地速 MAE 分别为 {row("C1","F").MAE_3_16:.4f}、{nt.MAE_3_16:.4f} m/s。移动管使名义恢复证明可用于这条参考轨迹，也在本例中带来可量化的跟踪代价。'))
    actkeys=[key for key in ac if key.endswith('_samples')]
    doc.append(p(f'全配置已保存网格中，检测到执行器幅值、速率或指令限幅的运行共 {int((ac[actkeys].max(axis=1)>0).sum())}/30 条，分项计数及积分持续时间见执行器表。所有保留状态的归一化 QP KKT 残差最大 {metrics.KKT_max.max():.2e}。提前退出组的位置 RMS 仅与共同网格上配对的 F 比较，静态组捕获从 0.05 s 开始。'))
    doc.extend(['<h2>论文可以据此回答的问题</h2>',p('R3-2 可以报告从原配平邻域到低速保持段的可恢复参考轨迹，并给出关键窗在实际飞行记录上的条件覆盖。C1 与 C2 相近的平台要求把适用范围限定为该预先设计的时间表及所测离散工况。C3 给出一条弱风记录，证书仍以时间窗和误差条件为前提。'),
    p('R3-4 可以给出相对于仅限推力的反演根域保护机制，并逐项量化收益来源。主要地速改善来自更换备份，同管端点增量为零。两次无管实验中，等预留攻角也避免了同一类上端根事件。端点力差的优势可以写为能用牛顿核算两侧预留与反馈损失，现有闭环结果支持两种保护在这些机动上效果接近。研究尺度尚无独立任务验收依据。原方案的较强判据未触发，至此完成确认阶段。'),
    '<h2>复算记录</h2>',p('第一轮 18 个 F/D/A 作业在积分前的风接口自检停止。自检将 14、20 s 的移动管与初始配平状态相配后调用了整个控制器。按原冻结风场直接核对真风、空速、反馈需求和端点力差后重新启动这 18 个作业；保留原自检、首轮日志与两版输入清单。全部 30 个积分使用相同的既定飞行配置和判据。独立 Python 复算核对惯性水平地速、气流三轴误差、原始反馈端点、解析风导和移动椭球成员。'),
    p('<a href="results/run_metrics.csv">运行指标</a>　<a href="results/contrasts.csv">双步长配对</a>　<a href="results/cross_case_F.csv">跨工况差</a>　<a href="results/common_windows.csv">共同时间窗</a>　<a href="results/events.csv">事件</a>　<a href="results/matched_replay.csv">同状态消融</a>　<a href="results/free_replay.csv">无管消融</a>　<a href="results/coverage.csv">证书覆盖</a>　<a href="results/actuators.csv">执行器</a>　<a href="README_复算.txt">复算说明</a>'),'</html>'])
    (H/REPORT).write_text('\n'.join(doc),encoding='utf8')

def package():
    files=sorted(p for p in H.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='release_manifest.json')
    hashes={p.relative_to(H).as_posix():digest(p) for p in files}
    manifest=H/'release_manifest.json';manifest.write_text(json.dumps(hashes,indent=2,ensure_ascii=False),encoding='utf8')
    pack=ROOT/'Aerospace_有限时长减速机动_C1-C3闭环确认_v1_0.zip'
    with zipfile.ZipFile(pack,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in [*files,manifest]:z.write(p,p.relative_to(H).as_posix())
    import hashlib
    with zipfile.ZipFile(pack) as z:
        assert z.testzip() is None and all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in hashes.items())
    print('Archive verified:',len(files)+1,'files;',pack.stat().st_size,'bytes',flush=True)
if __name__=='__main__':
    assemble();package()
