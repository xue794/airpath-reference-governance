"""Assemble the D0 development report from the paired saved runs."""
import html,zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from evaluate import *
def table(head,rows):
    def tr(tag,row):return '<tr>'+''.join(f'<{tag}>{html.escape(str(v))}</{tag}>' for v in row)+'</tr>'
    return '<table><thead>'+tr('th',head)+'</thead><tbody>'+''.join(tr('td',r) for r in rows)+'</tbody></table>'
def fig(tr):
    zh=FontProperties(fname='C:/Windows/Fonts/msyh.ttc',size=10)
    fig,ax=plt.subplots(2,1,figsize=(8.8,6.4),sharex=True)
    for a,col,style,txt in [('O','#65727b','-.','O 原静态备份'),('F','#11685d','-','F / D / A 移动备份（三者重合）'),('T','#ad4b29',':','T 无管，仅限推力等常规条件'),('N10','#384ea1','-','N10 无管，端点预留'),('Am','#9b671c','--','Am 无管，等预留攻角')]:
        x=tr[a,'005'];ax[0].plot(x.time,x.vground,style,color=col,label=txt,lw=1.65);ax[1].plot(x.time,x.fbHi,style,color=col,lw=1.65)
    x=tr['F','005'];ax[0].plot(x.time,x.vrequest,'k--',label='请求地速',lw=1);ax[0].axvspan(3,16,color='#e9eeee',zorder=-3);ax[1].axvspan(14,20,color='#e9eeee',zorder=-3)
    x=tr['T','005'].iloc[-1];ax[1].scatter([x.time],[x.fbHi],color='#ad4b29',s=18,zorder=8);ax[1].axhline(0,ls=':',lw=1,color='#506065')
    ax[0].set_ylabel('实际水平地速 (m/s)',fontproperties=zh);ax[1].set_ylabel('反馈上端力差 (N)',fontproperties=zh);ax[1].set_xlabel('时间 (s)',fontproperties=zh)
    for a in ax:
        a.spines[['top','right']].set_visible(False);a.grid(axis='y',alpha=.16);a.set_xlim(0,35)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',ncol=3,prop=zh,frameon=False)
    fig.tight_layout(rect=[0,0,1,.89]);fig.savefig(OUT/'D0_trajectories.png',dpi=160);plt.close(fig)
def main():
    result,contrasts,tr=run();fine=result[result.step=='005'].set_index('arm');fig(tr)
    abl=[];prefix=pd.read_csv(OUT/'prefix.csv')
    for s in ['005','010']:
        src=H/'runs';new=pd.read_csv(src/('N10_'+s)/'free_replay.csv')
        mask=new.NT_delta>1e-7;assert new[['T_KKT','Am_KKT']].max().max()<1e-8
        assert (new.loc[mask,'removed_phi15']<0).all()
        abl.append(dict(step=s,changed_against_T=int(mask.sum()),first_s=float(new.time[mask].min()),changed_against_Am=int((new.NA_delta>1e-7).sum()),max_input_difference=float(new.NT_delta.max()),min_upper_row=float(new.removed_phi15.min())))
        for m in ['D','A']:
            f=loadmat(src/('F_'+s)/'F.mat',simplify_cells=True);d=loadmat(src/(m+'_'+s)/(m+'.mat'),simplify_cells=True)
            assert np.array_equal(f['Xs'],d['Xs']) and np.array_equal(f['Ts'],d['Ts'])
    event=loadmat(H/'runs/T_005/T.mat',simplify_cells=True)
    old=ROOT/'contribution_probe_20261106/feedback/V31_C60_warm_fine_T_transfer.csv'
    ideal=float(pd.read_csv(old).Tstar_N.iloc[-1]);Tref=float(event['D'][-1,9]);Traw=float(event['D'][-1,10]);Tact=float(event['D'][-1,11])
    outcome=dict(runs=14,unique_arms=7,finished=int((result.cause=='completed').sum()),scope='D0 only',allow_confirmation=True,
        direct_speed_effect=bool(contrasts[0]['meets_scale']),same_tube_endpoint_effect=False,backup_free_endpoint_effect_vs_T=True,
        independently_resolved_no_tube_N10_vs_Am_N=float(contrasts[-1]['fine']),no_tube_ablation=abl,
        max_abs_prefix36=float(prefix.max_abs_36state_difference.max()),thrust_field_correction=dict(Tref=Tref,Tstar=ideal,Traw=Traw,Tact=Tact),
        confirmation_runs=0)
    assert outcome['finished']==10
    put('decision.json',outcome)
    main_rows=[]
    labels=['静态备份＋端点','移动备份＋端点','同管删端点','同管等预留攻角','无管推力对照','无管端点','无管等预留攻角']
    def fmt(x,p=3):return '—' if pd.isna(x) else f'{x:.{p}f}'
    for m,l in zip(ARMS,labels):
        rr=fine.loc[m];status={'completed':'完成','raw_command_speed_chart':'原始指令坐标域退出','feedback_upper_root':'反馈上端触界'}[rr.cause]
        main_rows.append([m+'　'+l,f'{rr.end_time:.6f}　'+status,fmt(rr.MAE_3_16),fmt(rr.vground_25_35),fmt(rr.position_RMS35),fmt(rr.min_fb_hi_window,4)])
    css='body{font:16px/1.66 "Source Sans 3","Microsoft YaHei",sans-serif;color:#173139;max-width:1040px;margin:32px auto;padding:0 25px}h1{font-size:29px}h2{font-size:21px;margin:28px 0 6px}p{margin:8px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:7px;text-align:left;border-bottom:1px solid #cdd6d7}th{background:#edf3f2}img{max-width:100%}a{color:#0c625c}.lead{border-left:4px solid #187068;background:#eef5f3;padding:8px 14px}small{color:#4e6268}'
    doc=['<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><title>D0 · 直线减速开发验证</title><style>'+css+'</style><main><h1>D0 直线减速的闭环验证与收益归因</h1><small>2026 年 11 月 8 日　｜　有限时长减速机动研究　｜　开发工况</small>',
      '<p class="lead"><strong>按预定规则，可以进入 C1–C3。</strong>换用已认证的移动备份后，减速窗的实际地速误差减少 1.21376 m/s，低速保持段平均达到 31.08973 m/s。在同一移动管内，完整端点、删行和等预留攻角三组给出相同的闭环轨迹。端点保护相对推力限制的作用出现在取消备份行的对照中。</p>',
      '<h2>七个配置与可比时间窗</h2><p>从 38.1 m/s 配平起飞行，3–11 s 按原五次函数将地速请求减至 31 m/s，随后保持至 35 s。航向固定、无风无噪，六自由度模型、下游控制器、QP 选择规则一致。3–16 s 地速指标为实际水平地速与请求地速绝对差的时间平均。表内为最大步长 0.005 s 的结果；另一路使用 0.01 s，各自记录间隔为 0.05 s。</p>',
      table(['配置','结束 / s','3–16 s 地速误差 / (m/s)','25–35 s 平均地速 / (m/s)','全程位置 RMS / m','14–20 s 最小反馈上端力差 / N'],main_rows),
      '<p>“—”表示没有覆盖整个指定窗口。位置 RMS 按原稿的规则在完整保存网格上计算；提前终止的轨迹另做共同保存时域比较。O 与 F 共有的 0.05–17.05 s 窗内，位置 RMS 分别为 28.8093 m 和 22.6755 m。无管 N10 的全程 RMS 为 17.6097 m，带管 F 为 38.3479 m，反映本机动的一项跟踪代价。当前无独立任务验收阈值。</p>',
      '<img src="results/D0_trajectories.png" alt="七种配置的实际地速与反馈上端力差"><p><small>灰带分别标出减速评价窗与反馈核查窗。横轴均为飞行时间，提前退出只画至最后有效记录。N10 与 Am 的曲线十分接近。</small></p>',
      '<h2>原静态备份与移动备份</h2><p>O 和 F 的主指标分别为 4.46507 与 3.25131 m/s。改善量在两个积分步长下均为 1.21376 m/s，超过预定的 0.606 m/s 研究尺度和 0.0005 m/s 数值漂移判据。O 在 17.088295 s 的首个边界为原始指令速度模降到 5 m/s；水平指令速度仍约为 4.96995 m/s，高于 2 m/s 的水平坐标门槛。此时四个反演端点力差均为正，QP 已有可行解。该边界描述制导指令坐标的适用范围。</p>',
      '<p>F 于 7.20 s 首次记录到原备份椭球外的参考，全程 557/701 个保存点在原域外，701/701 个点在新管内。该新管在 25–35 s 的整个截面与旧域分离；本次闭环实测的该段平均地速为 31.08973 m/s，治理参考换算的平均地速为 31.05949 m/s。原规划的 35–55 s 恢复尾段属于第 3 节存在性证明，本开发工况在 35 s 结束。</p>',
      '<h2>端点行与等预留攻角保护</h2><p>两种步长下，F、D、A 各自独立积分得到的 36 维状态在全部 701 个保存点逐值相同。在 F 的同一状态与同一名义输入下重解，删端点行和改用攻角界各有 0/701 个可辨指令变化，最大差均为零。F 的 614 个干预点均由备份行获得正对偶支持，端点行获得支持的点数为零。三组在关键窗内的最小上端反馈力差均为 23.27189 N，端点相对两种对照的闭环增量均为 0 N。</p>',
      '<p>去除备份行后，在 N10 的 701 个状态上改用 T，有 246 个点发生可辨输入变化，始于 13.85 s。这 246 个删行解均违反原上端点导数条件。T 在 14.652183 s 触及反馈上端，N10 和 Am 则均完成 35 s。14.65 s 的最后共同网格点上，T、N10、Am 的反馈上端余量分别约为 0.0722、13.3884、13.1231 N。N10 相对 Am 的关键窗最小裕度增量为 0.001973 N，通过步长漂移核查，远低于原定 1 N 量级。二者在 D0 上提供了接近的保护。</p>',
      '<h2>反馈盒覆盖与执行器</h2><p>在 F 的 14–20 s 关键窗，121/121 个保存点同时满足新管成员、三个分轴误差限值及 θF+θR≤0.70°。该窗姿态误差和最大约 0.02344°、所见上端力损失最大约 4.59029 N。这给出了第 3 节条件证书在本次飞行记录上的采样覆盖。实际误差是否连续驻盒，仍须连续时间误差界支持。</p>',
      '<p>F 的全程实际攻角为 4.8970–8.2966°，原始反馈推力为 7.0994–24.0701 N，实际推力为 7.1009–24.0673 N。七组在保留网格均未见所查执行器幅值或速率触限。其余各臂和各步长的力差、执行器、QP 残差与采样占用见附表。</p>',
      '<h2>复核与既有数字更正</h2><p>三组无管对照与 H31 对应历史轨迹按相同步长核对。在转向前的共同输入段或提前事件前，36 维状态最大绝对差小于 3.30×10⁻¹¹。额外用一个从零开始的短时风输入夹具检验真实飞机、风估计、反馈与事件读到同一风；D0 自身的真风与噪声始终为零，未取用 C3 的输出。每个保存点另行重建地速、三轴误差、原始反馈端点力差，并用有理数回算移动椭球成员。</p>',
      '<p>O 主积分在坐标守卫处停止，原始保存止于 17.05 s。自该点分别以 0.001/0.0005 s 最大步长回放，只为定位原 5 m/s 坐标边界放宽内部试算门槛，事件后不纳入任何比较。T 的上端触界由原四个事件函数直接定位。</p>',
      f'<p><b>第 2 节报告的推力字段更正。</b>H31/T 在 14.652183 s 的 {ideal:.6f} N 是 Tstar（理想反馈推力）。原始 MAT 文件的参考前馈推力 Tref 为 {Tref:.6f} N，原始反馈请求 Traw 为 {Traw:.6f} N，实际执行推力为 {Tact:.6f} N。此前文字把 Tstar 写成了参考推力。更正后，推力合规而反馈上端力差到零的判断保持成立。</p>',
      '<h2>开发判定与论文口径</h2><p>D0 满足进入确认工况的条件：无管对照复现端点作用，更换备份后地速跟踪也达到预定改善尺度。可以继续原定 C1、C2、C3，检验时长、低速目标和弱风变化。D0 支持 R3-2 所问的“低速域能否实际用上”；对 R3-4，它把已观察到的减速改善归于备份，给出了相对仅限推力的反演保护，也显示同管端点行与等预留攻角比较没有达到独立增益尺度。</p>',
      '<p><a href="results/run_metrics.csv">运行指标</a>　<a href="results/contrasts.csv">双步长配对差</a>　<a href="results/matched_replay.csv">同状态重解</a>　<a href="results/coverage.csv">证书覆盖</a>　<a href="results/events.csv">首次事件</a>　<a href="README_复算.txt">复算说明</a></p></main></html>']
    (H/'D0_直线减速闭环开发报告.html').write_text('\n'.join(doc),encoding='utf8')
    paths=sorted(p for p in H.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='release_manifest.json')
    manifest={p.relative_to(H).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    fp=H/'release_manifest.json';fp.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf8')
    pack=ROOT/'Aerospace_有限时长减速机动_D0闭环开发验证_v1_0.zip'
    with zipfile.ZipFile(pack,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in [*paths,fp]:z.write(p,p.relative_to(H).as_posix())
    with zipfile.ZipFile(pack) as z:
        assert z.testzip() is None and all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in manifest.items())
    print('Verified zip bytes',pack.stat().st_size)
if __name__=='__main__':main()
