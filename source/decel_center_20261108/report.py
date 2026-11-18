"""Build the Chinese report and verify the complete study archive."""
import html, json, zipfile, hashlib
from pathlib import Path
from html.parser import HTMLParser
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from evaluate_direct import H,PREV,OUT,sha,run

def table(head,rows):
    def tr(tag,row):return '<tr>'+''.join(f'<{tag}>{html.escape(str(x))}</{tag}>' for x in row)+'</tr>'
    return '<table><thead>'+tr('th',head)+'</thead><tbody>'+''.join(tr('td',row) for row in rows)+'</tbody></table>'
def fig(tr):
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.is_file():
        font_manager.fontManager.addfont(font);plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams['axes.unicode_minus']=False
    ff,axes=plt.subplots(2,2,figsize=(10.4,6.9),layout='constrained')
    for ax,sc in zip(axes[0],['C1','C2']):
        d=tr[sc,'005','Direct'];f=tr[sc,'005','F']
        o=pd.read_csv(PREV/'runs'/f'{sc}_O_005/replay.csv')
        ax.axvspan(3,16,color='#e9f0f2')
        for t,y,label,col,style in [(f.time,f.vrequest,'请求','#39444d','--'),(o.time,o.vground,'原静态 O','#996446',':'),(d.time,d.vground,'直接中心','#357c99','-'),(f.time,f.vground,'完整 F','#ba4729','-')]:
            ax.plot(t,y,style,color=col,lw=1.7,label=label)
        ax.set(title=sc+'   实际水平地速',xlim=(0,35),ylim=(30.5,38.7),ylabel='地速 (m/s)');ax.legend(frameon=False,ncol=2,fontsize=8)
    ax=axes[1,0]
    for sc,style in [('C1','-'),('C2','--')]:
        for arm,col,txt in [('Direct','#357c99','直接中心'),('F','#ba4729','完整 F')]:
            x=tr[sc,'005',arm];ax.plot(x.time,x.positionError,style,color=col,lw=1.6,label=sc+' '+txt)
    ax.set(title='实际位置到原请求的距离',ylabel='三维位置误差 (m)',xlim=(0,35));ax.legend(frameon=False,fontsize=8)
    ax=axes[1,1];sc='C1'
    ax.axvspan(14,20,color='#e9f0f2')
    for arm,col,txt in [('Direct','#357c99','直接中心'),('F','#ba4729','完整 F')]:
        x=tr[sc,'005',arm];ax.plot(x.time,x.fbHi,color=col,lw=1.9,label=txt)
        k=x[np.isclose(x.time,20)].iloc[0];ax.scatter([20],[k.fbHi],s=24,color=col)
    ax.set(title='C1   反馈反演上端力裕度',ylabel='上端力差 (N)',xlim=(14,26),ylim=(0,85));ax.legend(frameon=False,fontsize=9)
    for ax in axes.flat:
        ax.set_xlabel('时间 (s)');ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.18)
    ff.savefig(OUT/'center_comparison.png',dpi=155);plt.close(ff)
def build():
    tab,tr=run();fig(tr)
    cases=['C1','C2'];fine=tab.xs('005',level='step')
    cov=pd.read_csv(OUT/'coverage.csv',dtype={'step':str})
    dr=pd.read_csv(OUT/'paired_differences.csv',dtype={'step':str})
    sh=pd.read_csv(OUT/'speed_comparison.csv',dtype={'step':str}).query("step=='005'").set_index('case')
    dq=pd.read_csv(OUT/'cross_request_checks.csv')
    assert (tab.status=='completed').all() and (tab.samples==701).all() and pd.read_csv(OUT/'events.csv').empty
    assert (dq.max_abs_13plant_state_diff==0).all() and (dq.max_abs_4actuator_state_diff==0).all()
    assert ((cov.query("window=='critical'").J==121)).all()
    assert (tab[[c for c in tab if c.endswith('_samples')]]==0).all().all()
    def get(sc,a):return fine.loc[sc,a]
    def val(sc,metric):return dr[(dr['case']==sc)&(dr.metric==metric)].iloc[0]
    p1,p2=[sh.loc[sc,'fraction_of_O_F_gain_at_Direct']*100 for sc in cases]
    a,b=get('C1','Direct'),get('C1','F')
    body=['''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>直接跟随中心曲线 · C1、C2 六自由度闭环对照</title><style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#152e38;max-width:1020px;margin:2.5em auto;padding:0 1em}h1{font-size:1.85em}h2{margin:1.45em 0 .4em;font-size:1.22em}p{margin:.6em 0}table{display:block;overflow-x:auto;border-collapse:collapse;font-size:.91em}td,th{padding:.35em .65em;white-space:nowrap;text-align:right;border-bottom:1px solid #cfdcde}th{background:#eaf1f2}td:first-child,th:first-child{text-align:left}img{max-width:100%}a{color:#076758}</style><h1>直接跟随中心曲线能实现多少减速改善？</h1>''',
    f'<p>在静风 C1 和 C2 中，直接把已验证中心送入下游控制器，已分别获得原静态备份 O 到完整方案 F 的实际地速误差改善的 {p1:.2f}%、{p2:.2f}%。F 进一步降低了实际地速误差和全程位置 RMS，同时其关键窗最小反馈上端力裕度较直接中心少约 3.23 N。</p>',
    '<h2>对照如何运行</h2><p>本补充实验在已知 C1–C3 确认结果后提出，并在首条新轨迹积分前登记。两个请求均从 38.1 m/s 配平开始，3 s 起用 10 s 分别减至 31、31.5 m/s。直接组依据相同的 22 s 五次中心曲线给出空速及其解析导数，航迹角和航向保持水平直线，25 s 起中心空速为 31.25 m/s。两个请求的外层制导及误差积分继续记录，实际飞机接收同一解析参考。F 在相同时间表的移动管内根据原请求在线选取参考。使用同一六自由度模型、下游反馈、滤波、分配、执行器、初值、静风及无噪声设置。各自以 0.01 s 和 0.005 s 最大步长独立计算，输出间隔 0.05 s。</p>',
    '<p>四条新轨迹和四条相配的历史 F 轨迹都到达 35 s。均未观测到反演根事件或原始指令坐标域退出。对直接组，后者是被旁路的外层制导量，只作监视。以下主表给细步长数值，附表保存粗步长和每项差值。</p>',
    '<h2>地速与位置</h2>']
    body.append(table(['工况','方案','3–16 s 实际地速 MAE (m/s)','参考地速 MAE (m/s)','全程位置 RMS (m)','峰值位置误差 (m)','平台实际 / 参考地速 (m/s)','实际平台偏差 (m/s)'],
        [[sc,'直接中心' if arm=='Direct' else '完整 F',f'{r.actual_MAE:.4f}',f'{r.reference_MAE:.4f}',f'{r.position_RMS:.3f}',f'{r.position_peak:.3f}',f'{r.plateau_actual:.5f} / {r.plateau_reference:.5f}',f'{r.plateau_bias:+.5f}'] for (sc,arm),r in fine.iterrows()]))
    body.append(table(['工况','原静态 O 的 MAE','O 到 F 总改善','O 到直接中心','直接中心到 F','直接中心已取得的比例'],[[sc,f'{r.O:.4f}',f'{r.O_minus_F:.4f}',f'{r.O_minus_Direct:.4f}',f'{r.Direct_minus_F:.4f}',f'{100*r.fraction_of_O_F_gain_at_Direct:.2f}%'] for sc,r in sh.iterrows()]))
    body.append('<p>速度单位均为 m/s。比例为 (MAE<sub>O</sub>−MAE<sub>直接</sub>)/(MAE<sub>O</sub>−MAE<sub>F</sub>)，只解释这三组在 3–16 s 同一指标上的差值。O 在两个工况中均覆盖该窗，随后分别约在 18.154 s 和 19.208 s 退出。直接中心与 F 的速度差在两个步长均大于 0.0005 m/s 数值门槛，低于原方案 0.606 m/s 的研究尺度；位置 RMS 改善作为描述性结果。</p>')
    body.append(f'<p>把目标提高 0.5 m/s 后，直接组平台实际地速的变化为 0；F 的平台实际和参考地速分别增加 {val("F","plateau_actual").fine:.8f}、{val("F","plateau_reference").fine:.8f} m/s。按相同步长比较，直接组 C1、C2 的 13 维飞机状态与四个执行器状态在全部 701 点逐值相同。两组位置误差不同，是因为两条请求位置不同。F 的平台地速更接近 C1 的 31 m/s 请求；C2 的 31.5 m/s 请求在平台期由直接中心跟得更近。整段位置和减速窗速度与终端保持是不同指标。</p>')
    body.append('<img src="results/center_comparison.png" alt="两种请求的原静态、直接中心和完整方案实际地速；位置误差及反馈上端裕度"><h2>反馈反演裕度及证书覆盖</h2>')
    body.append(table(['工况','方案','窗口上端最小 (N)','窗口下端最小 (N)','全程上端最小 (N) / 时间 (s)','全程下端最小 (N) / 时间 (s)','原始推力请求范围 (N)'],[[sc,'直接中心' if arm=='Direct' else '完整 F',f'{r.fbHi_critical:.4f}',f'{r.fbLo_critical:.4f}',f'{r.fbHi_observed:.4f} / {r.fbHi_observed_t:.2f}',f'{r.fbLo_observed:.4f} / {r.fbLo_observed_t:.2f}',f'{r.Traw_min:.3f}–{r.Traw_max:.3f}'] for (sc,arm),r in fine.iterrows()]))
    body.append(f'<p>14–20 s 的上端最小值均出现在 20 s，下端最小值在 14 s。F 的上端最小值低 {val("C1","fbHi_critical").fine:.3f} N（C1）和 {val("C2","fbHi_critical").fine:.3f} N（C2），两者超过既定 1 N 研究尺度。F 的同窗下端最小值较直接中心高约 0.275 N。各组四个端点在记录中均保持正裕度；已保存网格未见执行器幅值或速率触限。两种方法、两个工况、两个步长在关键窗均有 121/121 点符合移动管、分轴误差盒及姿态条件。该数是采样覆盖。</p>')
    body.append('<h2>对审稿意见的含义</h2><p>R3-2 可报告同一可恢复中心时间表本身已实现约 89% 的所测减速改善。完整方案在管内选取参考仍能改善早期跟踪，并降低这两例的全程位置 RMS。直接组按设计对终值请求无反馈，F 的平台也几乎没有随请求移动，对这一固定时间表的适用范围应据此表述。</p><p>R3-4 可用这个基线量化采用完整在线参考选择所获的跟踪增量及所用掉的上端裕度。比较同时改变了参考生成方式，端点行的独立贡献仍由先前同状态消融和无管根事件检验。89% 只适用于两个已见静风请求的同窗差值分解；四条补充轨迹不用于推断弱风、随机任务成功率或者其他治理器排名。</p>')
    body.append('<h2>复核与文件</h2><p>在配对 F 的 7 个时刻，用同一飞机状态和同一解析中心输入核对下游反馈、滤波、执行器与风估计输出，得到逐值一致。中心解析导数在 0–35 s 的分片和拼接点也与原移动备份实现吻合。Python 逐点从原始状态独立重建实际与参考地速、位置误差及两端反馈力差，并与历史 F 汇总核对。首次预检发现接口仍按原参数个数读默认配置，在积分前修正并重新锁定。初始登记、修正记录和全部原始轨迹均在包内。</p>')
    body.append('<p>'+ '　'.join(f'<a href="{url}">{name}</a>' for url,name in [('results/run_metrics.csv','逐运行指标'),('results/paired_differences.csv','双步长配对表'),('results/speed_comparison.csv','速度改善分解'),('results/coverage.csv','证书采样覆盖'),('results/events.csv','事件'),('results/cross_request_checks.csv','同中心状态核查'),('protocol.json','积分前登记'),('README_复算.txt','复算说明')])+'</p></html>')
    report=H/'直接跟随中心曲线_闭环对照报告.html';report.write_text('\n'.join(body),encoding='utf8')
    class Check(HTMLParser):
        def handle_starttag(self,tag,attrs):
            for k,v in attrs:
                if k in ['href','src']:assert (H/v).is_file(),v
    Check().feed(report.read_text(encoding='utf8'))
def package():
    files=sorted(p for p in H.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='release_manifest.json')
    hashes={p.relative_to(H).as_posix():sha(p) for p in files}
    manifest=H/'release_manifest.json';manifest.write_text(json.dumps(hashes,indent=2,ensure_ascii=False),encoding='utf8')
    pack=H.parent/'Aerospace_有限时长减速机动_直接中心闭环对照_v1_0.zip'
    with zipfile.ZipFile(pack,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in [*files,manifest]:z.write(p,p.relative_to(H).as_posix())
    with zipfile.ZipFile(pack) as z:
        assert z.testzip() is None and all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in hashes.items())
    print('Verified archive',pack.name,len(hashes)+1,pack.stat().st_size)
if __name__=='__main__':build();package()
