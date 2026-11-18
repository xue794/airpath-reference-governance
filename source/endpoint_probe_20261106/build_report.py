"""Summarize retained endpoint probes; never compare unequal-horizon RMS."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H=Path(__file__).resolve().parent
F=H/'full6dof'
reference=json.loads((H/'reference/continuous_summary.json').read_text())
raudit=json.loads((H/'reference/continuous_audit.json').read_text())
static=json.loads((H/'static/static_summary.json').read_text())
validation=json.loads((H/'review/final_validation.json').read_text())
assert validation['all_expected_results_present'], 'Final validation is incomplete'
assert validation['all_available_run_checks_pass'], 'Final validation checks failed'
groups={}
for speed in (31,32):
    p=F/f'V{speed}_C60_warm_fine'
    assert (p/'all_summary.json').exists(), f'Pending full-case results: {p}'
    groups[speed]={r['method']:r for r in json.loads((p/'all_summary.json').read_text())}

names={'T':'Thrust only','N0':'+ endpoints, 0 N reserve','N10':'+ endpoints, 10 N reserve'}
colors={'T':'#d55e00','N0':'#0072b2','N10':'#009e73'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.grid':True,'grid.alpha':.18,'savefig.facecolor':'white'})

fig,axs=plt.subplots(2,2,figsize=(12.8,8),constrained_layout=True)
for method,rname in [('T','thrust15'),('N0','full17_N0'),('N10','full17_N10')]:
    a=np.genfromtxt(H/f'reference/continuous_{rname}_maxstep0.005.csv',delimiter=',',names=True)
    axs[0,0].plot(a['t'],a['rawNhi'],color=colors[method],lw=1.7,label=names[method])
    axs[1,0].plot(a['t'],a['Tff'],color=colors[method],lw=1.7)
    if method=='T':
        axs[0,0].plot(a['t'][-1],a['rawNhi'][-1],'x',color=colors[method],ms=8,mew=2)
        axs[1,0].plot(a['t'][-1],a['Tff'][-1],'x',color=colors[method],ms=8,mew=2)
    b=np.genfromtxt(F/f'V31_C60_warm_fine/{method}_samples.csv',delimiter=',',names=True)
    axs[0,1].plot(b['time_s'],b['feedback_raw_upper_N'],color=colors[method],lw=1.7)
    axs[1,1].plot(b['time_s'],b['position_error_m'],color=colors[method],lw=1.7)
    if groups[31][method]['status']!='completed':
        axs[0,1].plot(b['time_s'][-1],b['feedback_raw_upper_N'][-1],'x',color=colors[method],ms=8,mew=2)
        axs[1,1].plot(b['time_s'][-1],b['position_error_m'][-1],'x',color=colors[method],ms=8,mew=2)
titles=['A  Reference-only: upper inversion margin','B  Full 6DOF: feedback upper inversion margin',
        'C  Reference-only: thrust already constrained','D  Full 6DOF: position error (truncated at domain exit)']
for ax,title in zip(axs.flat,titles):
    ax.set(title=title,xlabel='Time (s)',xlim=(0,35))
    ax.axvspan(16,19,color='gray',alpha=.07,lw=0)
axs[0,0].set_ylabel('N');axs[0,1].set_ylabel('N');axs[1,0].set_ylabel('N');axs[1,1].set_ylabel('m')
for ax in axs[0,:]:ax.axhline(0,color='#333333',ls='--',lw=1)
axs[0,0].axhline(10,color='#777777',ls=':',lw=1)
axs[1,0].axhline(3,color='#333333',ls='--',lw=1);axs[1,0].axhline(57,color='#333333',ls='--',lw=1)
axs[1,0].set_ylim(0,60)
axs[0,0].legend(fontsize=9,loc='upper right')
fig.suptitle('Incremental endpoint protection with the SAME coupled thrust constraints\n'
             '31 m/s request; no backup-set row; crosses denote terminal inversion-domain events',fontsize=13)
fig.savefig(H/'endpoint_evidence.png',dpi=180);fig.savefig(H/'endpoint_evidence.svg');plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(12.8,7.6),constrained_layout=True)
for method in ('T','N0','N10'):
    b=np.genfromtxt(F/f'V32_C60_warm_fine/{method}_samples.csv',delimiter=',',names=True)
    for ax,field in zip(axs.flat,['ref_raw_upper_N','feedback_raw_upper_N','Traw_command_N','position_error_m']):
        ax.plot(b['time_s'],b[field],color=colors[method],lw=1.7,label=names[method])
        if groups[32][method]['status']!='completed':ax.plot(b['time_s'][-1],b[field][-1],'x',color=colors[method],ms=8,mew=2)
for ax,title,unit in zip(axs.flat,['Reference upper inversion margin','Feedback upper inversion margin','Raw thrust command','Position error'],['N','N','N','m']):
    ax.set(title=title,xlabel='Time (s)',ylabel=unit,xlim=(0,35));ax.axvspan(16,19,color='gray',alpha=.07,lw=0)
for ax in axs[0,:]:ax.axhline(0,color='#333333',ls='--',lw=1)
axs[0,0].legend(fontsize=9)
fig.suptitle('Complementary full-6DOF check: 32 m/s request\nSame controller and thrust limits; terminal-domain events are not aircraft crashes',fontsize=13)
fig.savefig(H/'full6dof_32mps.png',dpi=180);plt.close(fig)

fulltable=[]
for speed,g in groups.items():
    for m in ('T','N0','N10'):
        a=g[m];event=a['event_indices'];event=[event] if isinstance(event,(float,int)) else event
        status='完成35秒' if a['status']=='completed' else '到达'+('/'.join({1:'参考下端点',2:'参考上端点',3:'反馈下端点',4:'反馈上端点'}.get(int(i),'未知事件') for i in event) or a['status'])
        fulltable.append(f"{speed} m/s | {m:3} | {a['end_time_s']:.9f} s | {status} | min ref upper {a['min_ref_raw_upper_N']:.6f} N | min feedback upper {a['min_feedback_raw_upper_N']:.6f} N")
ref=[r for r in reference if r['maxStep']==.005]
ref_table='\n'.join(f"{r['method']} | end {r['endTime']:.9f} s | Tff [{r['minT']:.6f},{r['maxT']:.6f}] N | min raw N+ {r['minRawNhi']:.9f} N" for r in ref)
g31=groups[31]
text=f'''已有推力约束之外，端点/反演域维护还能增加什么收益？
第二轮探索性验证报告，2026-11-06

结论
已经找到端点约束相对于“同一耦合推力受限基线”的独立增量作用：
推力保持3–57N时，参考更新仍可能使所选0–9度反演区间失效；显式
端点约束能够在本例中阻止这种参考域退出。零端点预留也有这个作用，
所以收益不完全来自人为10N预留。

完整6DOF层的结论更有限：参考域维护不能自动保证反馈指令始终留在
反演区间内。31m/s例中T和N0在减速阶段先到反馈上端点，N10完成35秒；
N10仍有较大的位置误差。本轮是接口可实现性的有限证据，不是物理
全局安全、优于所有传统控制方法或满足任务精度要求的证明。

一、基线及范围
T = 原仿射核的[1:13,16:17]行：速度/爬升/率/正法向/银行角，另加
相同耦合Tff的3–57N推力约束。T已经使用耦合隐式根与隐式导数。
N0 = 原1:17行，将DeltaN设为0，只维护原反演区间。
N10 = 原1:17行，保留投稿中的DeltaN=10N额外预留。
所有方法使用相同QP、尺度、推力预留、控制器、执行器、初值和请求。
各组均不施加KB，目的是排除小备份集遮蔽效应；都没有继承原定理的
递归可行保证。这些新结果也不等于原18行Scalar方案在其认证域内的
增益已经被验证。

二、静态证据
固定450个状态格点，277个同时满足常规状态、psi1、反演存在、推力
范围和10N端点预留。每点7个请求，共1939次同状态QP比较，全部成功。
T有629次违反10N端点导数条件；583次连零预留端点导数条件也违反。
这些是有限探索格的计数，不是随机总体概率。
简单例：V=31m/s，水平、各参考率为0，alpha=8.094339度，Tff=24.782992N。
请求gamma=8度，T给出的u_gamma=0.2004054rad/s²，原端点导数条件
为-120.542149N/s；N0/N10分别限制至0.0285879/0.0143262rad/s²。
当前根存在且推力合法，并不能使端点行成为冗余约束。

三、从合法配平出发的连续参考动力学
请求：3–11秒C2减速至31m/s；16–19秒C2转向60度；总35秒。
保留12个有限筛查请求的全部结果（8个命中、4个阴性）；筛查使用
采样保持，仅作找例，不用于连续安全结论。选择公开披露为探索性。
后续连续ODE分别采用MaxStep=.01/.005，RelTol=1e-9、AbsTol=1e-11。

{ref_table}

T在17.522554707秒到达原始上端点，此时Tff=29.902329726N，bank约
20.1307度，course rate约6.78354度/秒。端点裕度导数=-49.780860324N/s，
明确向外，不是只与边界相切。所有保存点的推力预留、常规状态及
psi1条件均通过数值检查。两步长的事件时间差为5.77e-10秒。
N0/N10都完成35秒，最低原上端点裕度约0.00923346/10.00001519N。
沿样本原反演区间的映射导数系数下界为正；该检查使用普通浮点，
不是新的外向舍入区间认证。

四、完整六自由度结果（细步长）
采用原W2-D的无风、同名义飞机模型、无测量噪声配置；本轮未检验
风扰动、模型失配或传感器噪声下的鲁棒性。
以完全相同的C2速度/航向构造地面轨迹，位置由同一速度律积分。
实际参考指令仍来自原位置制导和估计器，不能把参考层与6DOF的同一
名义请求误认为具有相同内部轨迹。

{chr(10).join(fulltable)}

31m/s中T/N0的边界事件均早于16秒转向开始，是减速阶段的反馈域退出。
参考根仍存在、参考推力仍在预留范围内；这直接说明局部裕度传递
不能被省略。N10完成35秒，但全时域采样位置RMS为
{g31['N10']['RMS_until_stop_m']:.6f}m，不能自动称作可接受的工程跟踪。
32m/s相邻请求得到相同类型的结果：T/N0仍先到反馈上端点，N10完成
35秒，全程位置RMS为{groups[32]['N10']['RMS_until_stop_m']:.6f}m。该组也没有消除跟踪误差偏大的问题。
T/N0的截断RMS绝不与N10全程RMS直接排名。共同区间积分RMS、两步长
逐状态差异、所有接受点的域检查见review/final_validation.json/.txt。
31m/s的共同区间0–14.652183秒，三组RMS均约9.240m；目前额外收益是
保持接口可实现性，本试验并未显示明显的共同区间跟踪精度优势。
两组请求、三种配置、两种MaxStep共12条完整核查记录齐全，保存点的
原始根裕度、事件截断、常规参考约束和推力盒检查均通过数值容差。

五、数值事件与计算核查
根延拓只为ODE试探点定位终止事件：先调用原0–9度根，失败时用更宽
区间临时延拓，端点条件始终按原0–9度计算。参考/反馈上下端点共4个
terminal事件；有效轨迹止于首次原域边界。没有将延拓后的飞行当结果。
原始粗步长使用原QP数值起点；加速版只改变FeasibleStart提议，仍由
同一历史w2b_qp3核检查可行性和KKT。独立同轨迹比较与步长核查已保留。
原冗余细化批次中尚未完成的重复N10被停止，全部科学case继续由相同
QP的加速版完成。不是删除负结果，记录在full6dof/REDUNDANT_BATCH_STOPPED.txt。
原V31 coarse JSON的failed_trial_time_s曾误标最后RHS时间，原文件保留，
事件时间一律用event_times_s；后续文件已区分这两个量。

六、不能越过的贡献边界
1. 此处验证的是端点/定义域维护相对于推力限制的独立作用。T基线本来
   就使用耦合模型，不能写成“证明耦合反演优于解耦反演”。
2. 在M_alpha>0且DeltaN=0时，端点状态集等价于alpha根在0–9度内。
   在攻角上界CBF增益同取1/s的静态检查中，标准攻角CBF也会拒绝T的
   更新。内部导数限制因状态相关缩放
   不必完全相同，但没有证据表明本方案优于“推力＋攻角CBF”方法。
3. 9度是所选控制反演区间的上界，代码更宽攻角研究域上界为10度。
   到达9度域边界不是飞机已失速、坠毁或违反完整物理飞行包线。
4. 新轨迹远离原认证KB。要纳入正文作为方法的可用范围证据，还须
   解决备份集覆盖与反馈裕度传递，不能直接沿用原18行定理。
5. 没有全局安全保证、任务精度阈值或统计泛化声明。采样无违例不证明
   采样间全部物理状态安全。所有停止、阴性与数值舍入结果均保留。

七、对论文修改的直接意义
可以保留一个具体、可复现的主张：仅满足推力限制并不足以维持所用
耦合反演的定义域；端点约束能在推力远离上下限时发挥独立作用。
建议把这一点作为专用约束构造的机制证据，配合其实现代价和适用边界。
新颖性方面仍需与带攻角/载荷保护的强基线逐项比较。工程方面需要
让可验证备份构造覆盖有收益的区域，并定量连接参考裕度与反馈实现。

所有原始代码、第一轮归档和历史结果保持不变。
复现命令见REPRODUCE.txt；全部搜索与轨迹数据均在本目录。
'''
(H/'RESULTS_CN.txt').write_text(text,encoding='utf-8')

# Confirm the previously delivered first-round artifacts have not changed.
first=H.parent/'revision_probe_20261106'
items=json.loads((first/'artifact_manifest.json').read_text(encoding='utf-8'))
changed=[]
for item in items:
    p=first/item['path']
    if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:changed.append(str(p))
assert not changed, f'First-round artifacts changed: {changed}'
(H/'first_round_integrity.json').write_text(json.dumps({'files_checked':len(items),'all_unchanged':True},indent=2))
print(json.dumps({'report':str(H/'RESULTS_CN.txt'),'first_round_unchanged':len(items),'fine_full_cases':6}))
