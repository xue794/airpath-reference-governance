from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
summaries=json.loads((root/'continuous_summary.json').read_text())
audit=json.loads((root/'continuous_audit.json').read_text())
screen=json.loads((root/'screening_summary.json').read_text())
fine=[s for s in summaries if s['maxStep']==.005]
colors=['#d55e00','#0072b2','#009e73']
names=['Thrust only (15 rows)','+ endpoint, no reserve','+ endpoint, 10 N reserve']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
fig,axs=plt.subplots(2,2,figsize=(12.4,8),constrained_layout=True)
for q,color,name in zip(fine,colors,names):
 a=np.genfromtxt(root/f"continuous_{q['method']}_maxstep0.005.csv",delimiter=',',names=True)
 for ax,field in zip(axs.flat,['rawNhi','Tff','V','chi']):
  y=a[field];y=y*180/np.pi if field=='chi' else y
  ax.plot(a['t'],y,label=name,color=color,lw=1.8)
axs[0,0].axhline(0,color='black',ls='--',lw=1);axs[0,0].axhline(10,color='gray',ls=':',lw=1)
for bound in (3,57):axs[0,1].axhline(bound,color='gray',ls='--',lw=1)
for ax,title,ylabel in zip(axs.flat,['Original upper endpoint margin','Coupled reference thrust','Reference airspeed','Reference course'],['N','N','m/s','deg']):
 ax.set(title=title,xlabel='Time (s)',ylabel=ylabel,xlim=(0,35));ax.axvspan(16,19,color='gray',alpha=.08,lw=0)
axs[0,0].set_ylim(-5,155);axs[0,1].set_ylim(0,60);axs[1,1].set_ylim(-2,65)
axs[0,1].legend(loc='upper right',fontsize=9)
fig.suptitle('Reference-level endpoint ablation beyond identical thrust limits\nContinuous feedback; original-domain event stops thrust-only run',fontsize=13)
fig.savefig(root/'reference_endpoint_validation.png',dpi=180);fig.savefig(root/'reference_endpoint_validation.svg');plt.close(fig)
lines=['第二轮：相同耦合推力约束下，反演端点约束的增量作用','',
'范围：仅参考动力学，不是完整飞机物理安全验证。原第一轮代码和结果未修改。',
'基线thrust15使用1:13及16:17，保留相同耦合Tff、3–57 N预留、QP权重和求解器。',
'对照full17_N0添加原0–9deg反演区间端点条件；full17_N10再添加原稿10 N端点预留。',
'三组均无备份行，因此均不继承原备份集递归可行性证书。',
'基线已采用耦合模型；本轮隔离端点/定义域保护的增量，不比较耦合与非耦合模型。','',
'有限探索：速度[30,31,32,34]×course[30,60,90]共12个请求，全部从合法38.1m/s配平开始。',
'请求为C2：3–11秒减速，16–19秒改变course，35秒结束。screen步长0.025s。',
'8个筛查请求到达原端点，4个完成，全部结果保留，不能把比例解释为发生概率。',
'粗采样中速度30/31例出现轻微推力预留越界，因此不把筛查作为连续保证；另做连续验证。','',
'连续验证锁定速度31m/s、course60deg，ode45 RelTol=1e-9 AbsTol=1e-11，',
'MaxStep分别0.01/0.005s，输出采样0.02s，原raw endpoint=0为终止事件。',
'仅为根事件定位，对越过原区间的RHS试探点用[-2,12]deg根搜索光滑延拓；',
'原0/9deg端点约束c/G始终不变。延拓点不计为原模型的有效轨迹。','',
'最细步长结果：']
for s in fine:
 lines.append(f"{s['method']}: {s['status']}, t_end={s['endTime']:.9f}s, Tff=[{s['minT']:.9f},{s['maxT']:.9f}]N, rawNhi_min={s['minRawNhi']:.9g}N, thrust_violations={s['thrustReserveViolations']}, conventional_violations={s['conventionalViolations']}, extended_RHS_calls={s['extensionCalls']}")
lines+=['','端点事件与收敛审计：',json.dumps(audit,ensure_ascii=False,indent=2),'',
'映射导数下界使用整个原0–9deg区间的保守系数界，普通浮点求值，非外向舍入认证；不得套用于负alpha扩展区间。',
'意义：沿连续参考动力学，从合法初值确实可在推力预留仍成立时到达反演区间上端点。',
'此时原系统不能继续保持既定的严格内部反演定义域；负的端点裕度导数显示是向外运动而非擦边。',
'添加零预留端点约束即可在本轨迹保持原反演区间；10N版本另提供储备，因此收益不只是人为10N裕度。',
'这不是实际飞机失稳，不是对所有AoA/载荷受限常规包线方法的优势证明，也不是全局递归可行证明。',
'不同终止时域的RMS不得直接比较；CSV/JSON中的post16指标保留原终止时域，仅用于各轨迹描述。',
'阴性结果：速度32/course30与速度34的三种course请求均未触端点。','',
'复现：在MATLAB R2018a加本目录路径，依次运行run_endpoint_search、run_endpoint_continuous、audit_endpoint_validation；',
'之后用Python运行summarize_endpoint_reference.py。prepare_endpoint_kernel.py只从冻结源码生成唯一名称的事件延拓核。']
(root/'REFERENCE_RESULTS.txt').write_text('\n'.join(lines),encoding='utf-8')
print('Reference report and figure generated')
