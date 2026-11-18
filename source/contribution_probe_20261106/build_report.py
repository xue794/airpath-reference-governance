"""Assemble the final third-round scientific interpretation and provenance."""
from pathlib import Path
import hashlib,json

H=Path(__file__).resolve().parent
comparison=json.loads((H/'comparison_audit.json').read_text())
assert comparison['complete'] and comparison['all_available_run_checks_pass']
assert all(r['same_status'] and r['same_events'] for r in comparison['refinement'])
aoa=json.loads((H/'aoa/aoa_static_summary.json').read_text())
backup=json.loads((H/'backup/backup_summary.json').read_text())
transfer=json.loads((H/'feedback/transfer_summary.json').read_text())
posthoc=json.loads((H/'feedback/posthoc_attitude_summary.json').read_text())
assert len(transfer)==6
assert len(posthoc)==2
ft={int(x['case'][1:3]):x for x in transfer if x['method']=='N10'}
rows=[]
for p in comparison['matched_comparisons']:
    if p['methods'][0]!='Am':continue
    m=p['metrics']
    rows.append(f"{p['speed']} m/s | N10: {m['N10']['common_interval_RMS_m']:.6f} m | Am: {m['Am']['common_interval_RMS_m']:.6f} m | 两组均35秒")
trrows=[]
for v in (31,32):
    s=ft[v]
    trrows.append(f"{v} m/s | 最小epsilon={s['min_endpoint_radius_mps']:.6f} m/s | 逐点覆盖 {s['endpoint_tube_inclusion_samples']}/{s['samples']} ({100*s['endpoint_tube_inclusion_samples']/s['samples']:.2f}%)")
refine_pos=max(x['max_position_vector_difference_m'] for x in comparison['refinement'])
refine_time=max(x['end_time_difference_s'] for x in comparison['refinement'])
report=f'''专用约束能否成为论文核心贡献？第三轮小规模研究报告
2026-11-06

判断
现有证据尚不足以把“端点约束本身比攻角保护更强”作为论文核心贡献。
同耦合模型、同允许集合和同预留的攻角保护，复现了本轮端点保护的
主要闭环收益。端点形式有可计算性和求根前检查的实现特点，但本轮
没有证明它具备独立的性能或同等优化实现下的计算优势。

在保证方面，已经得到新的严格非零局部参考备份证书，并定量算出了
反馈裕度传递的条件范围；但这些结果尚未覆盖产生动态收益的整段任务，
更不能构成新的主动约束/饱和阶段闭环误差保证。
这不是证明研究不可能成立，而是定位了现稿与可成立贡献之间的缺口。

一、强攻角对照的定义与独立核验
N0/N10：第二轮原端点力裕度0/10N，复用已核查结果。
A0：原耦合隐式攻角根，固定0..9度上下界CBF，k=1/s。
Am：同一个根，但上下界随状态变化，满足
 M(alphaLower,V,X)=M(alphaLo,V,X)+10N，
 M(alphaUpper,V,X)=M(alphaHi,V,X)-10N。
边界的V/X导数完整进入CBF，包括加速度uV，不把边界冻结。
两条新攻角行替换原14:15行；原1:13和16:17、推力限值、QP目标/尺度、
原求解器、控制器、执行器和参考波形不变。三类比较均无KB行。
当前攻角根直接复用推力核的结果，没有让A0人为多求一次根。

在严格单调且正则的反演区间内，N0与A0描述同一状态集；N10与Am
也描述同一状态集。端点函数等于正的状态相关割线斜率乘攻角间隙；
因此在边界的向内条件相同，但内部固定增益CBF不必相同。
319个有根旧格点全部通过状态集合符号一致性核验；277个合法状态
上的3324项中心差分核验，最大解析导数误差约8.05e-11。
1939组同状态请求中，A0/Am QP全部成功。与N0/N10相比，最优更新
分别有585/641组不同；交叉违反导数行的情况两个方向都有。
这排除了“公式完全相同”和“一种形式处处包含另一种”的误述。

Am是本轮根据力裕度定义构造的强匹配诊断基线，并非一篇现有文献的
直接复现。重参数化关系本身不能证明已有研究已经提出同一构造，也
不能单凭这一点否定全部新颖性；它检验的是功能与性能上的独立增益。

二、完整六自由度闭环结果
无风、同名义模型、无测量噪声，35秒请求。31m/s比较A0/Am，32m/s
比较Am；各用MaxStep=.01/.005，新增六条运行，全部检查通过。
粗细终止类型一致，最大终止时间差{refine_time:.3g}秒，最大共有保存点
位置向量差{refine_pos:.3g}m。所有接受点止于首次原反演域事件。

31m/s：A0在15.374137862秒、N0在15.420596393秒到反馈反演上边界。
两者参考反演仍合法，不能把这个结果描述为飞机已经失速。
同一共同区间内积分位置RMS分别9.475252/9.475192m，几乎相同。

10N匹配比较（完整0–35秒、误差平方时间积分RMS）：
{chr(10).join(rows)}

这里使用统一的时间积分RMS；它与前轮直接对701个等间隔样本取均值的
RMS有细小差别。直接样本RMS也保存在各原始summary.json中。
本轮没有看到端点形式更好的跟踪表现；该有限样本也不意味着Am在所有
任务上更优。两种形式都有约42–47m的全程RMS，任务可接受性尚未建立。
状态集合匹配后基本复现收益，说明上一轮“可完成仿真”的改善主要由
域保护和预留解释，不能归因于端点表达形式的独有功能。

三、备份证书：能在局部建立，但未覆盖动态收益
复用原10^18定点整数外向区间核，验证了有限新候选的非零参考邻域。
在本有限搜索中，已通过的代表值为：
  31m/s、零转向率：rho=.003；
  32m/s、零转向率：rho=.03；
  34m/s、零转向率：rho=.05；
  34m/s、转向率中心.15rad/s：rho=.005。
它们保留10N端点及3–57N推力预留、初始化、正则性及操作CBF条件。
转向率中心改变后的区间核和相对平衡点反馈已单独列明；它不是实际
飞机恒转弯证书，也没有证明各局部集合之间能安全切换或恢复。

原KB(rho=.1)在两条完整N10轨迹中仅覆盖112/701和114/701个保存点。
上述新候选在这两条完整轨迹上的新增覆盖为0；所以并集覆盖仍约16%。
两条轨迹所需max VB分别20.2552和19.1301，远大于原rho=.1。
原固定备份满足全部1:17操作行的保存点只有183/701和195/701。
例如31m/s、16.80秒状态下，若采用原备份，N_plus导数行约-81.9836N/s；
该负号已由外向区间独立确认。这是在审计假设备份控制，并非保存
轨迹实际执行了这条不合法控制，也不是飞机物理失效。

参考层轨迹的候选并集覆盖改善到535/1751点，但关键17.52秒状态距
最近候选仍为VB/rho约15.97。31m/s静态机制例可以位于新的局部证书
中心，这说明有收益相关局部区域可验证，但不能替代动态任务覆盖。
还保留了严格内部反例：31m/s中心、rho=.01的椭球内有点不满足10N
端点预留；rho=.005仅在限定节点数内未完成区间覆盖，不可称为不可行。
未搜索或否定所有椭球形状、增益、轨迹备份及有限时域恢复构造。

四、反馈传递：从定性问题变成了可核查数量
直接从所用反馈反演推导两端点的精确差分：
 h_lower,c=h_lower,f+DeltaNforce-DeltaM_lower，
 h_upper,c=h_upper,f+DeltaM_upper-DeltaNforce。
这里DeltaNforce是反馈法向力需求与参考法向力需求之差，不能与固定
10N预留参数混淆。六条旧轨迹的精确差分/保存量重构误差均在约1e-12N内。

失败T/N0轨迹的关键误差是约-0.42至-0.44度的爬升角跟踪误差，它贡献
约6N上端点裕度损失；该阶段尚未开始转向。参考有根不自动保证反馈有根。

对每个保存参考点，采用明确的联合误差盒
 |eV|<=epsilon，|eGamma|,|eChi|<=epsilon/38.1（角度为rad）。
区间计算给出满足两端点条件的局部充分半径，取值搜索上限为1m/s；
到达上限只表示已验证到此，不能称真实最大半径。
N10结果如下：
{chr(10).join(trrows)}

88%/91%是逐点误差盒包含率，不是飞机安全概率，也不是连续时间覆盖率。
盒外的实际保存点仍可能安全；充分界未覆盖不能推断物理违例。
这些盒沿已知参考轨迹事后计算；尚无主动干预和饱和阶段的先验闭环
误差界保证误差始终留在盒内，因此不能写成新的闭环安全定理。

反馈原始推力命令Traw还依赖姿态实现误差引起的投影因子kappa，仅用
三个航迹误差不能独立保证它。反馈报告另列姿态条件及物理0–60N检查。
原0.70度仅作诊断；补充2.5度阈值由观察最大值向上取整，仍非先验保证。
在2.5度姿态条件下，两端点及0–60N原始推力指令的联合局部误差盒，
对31/32m/s分别覆盖390/701和409/701个保存点；这仍只是所选充分盒
的事后包含率。全部观察姿态误差小于2.5度不等于证明以后始终如此。
两条N10轨迹Traw最低约2.269/2.009N，低于参考预留3N但仍在物理0–60N
范围内。这是预留被使用，不是原物理推力盒违例，也不能拿它否定预留设计。
执行器限幅前的原始反馈推力指令与实际执行器推力始终分开记录。

五、对论文定位的建议
当前不宜将核心主张写成“提出了攻角保护无法实现的反演安全能力”，
或“在实用任务范围内已获得端到端递归可行与鲁棒安全”。
可以保留的具体价值是：端点形式能在求根前检查反演括号和力预留，
并以显式端点映射/偏导组织可计算约束；本轮也给出了局部证书与误差
传递的定量工具。这些是值得进一步发展的实现与验证贡献候选。

要把它提升为核心贡献，需要解决至少一个实质问题并形成完整证据：
1. 证明这种表达在同等优化的强攻角基线下，显著降低认证/实现成本，
   或能认证对照难以处理且实际有用的范围。本轮没有计时或认证成本
   公平比较，不能拿Am当前多求两个根就宣布所有攻角方法更慢。
2. 设计覆盖收益任务的状态相关/轨迹恢复备份，并验证过渡和递归可行；
   再把反馈误差界接到量化的反演/物理约束裕度，而非只检查参考层。
仅新增几组Monte Carlo、加引用或调整贡献措辞，无法填补这些缺口。
新的核心贡献是否具有文献新颖性仍需逐项对照原文；本轮是机制和
可验证性实验，不是完成了新的系统文献检索。

本轮材料说明
aoa/：强基线公式、静态/导数检查、源代码与反馈数学独立审读。
full6dof/：六条新增运行的MAT、CSV、JSON及冻结派生源。
backup/：严格证书、未完成候选、严格反例、保存点覆盖和完整脚本。
feedback/：精确差分、误差盒半径、姿态条件、失败归因和全部轨迹记录。
所有实验均为有限探索，保留负结果；完整复现见REPRODUCE.txt。
原两轮数据及归档保持不变，哈希核验见prior_round_integrity.json。
'''
(H/'RESULTS_CN.txt').write_text(report,encoding='utf-8')
integrity=[]
for name in ['revision_probe_20261106','endpoint_probe_20261106']:
    folder=H.parent/name;items=json.loads((folder/'artifact_manifest.json').read_text(encoding='utf-8'))
    changed=[]
    for item in items:
        p=folder/item['path']
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:changed.append(item['path'])
    assert not changed,(name,changed)
    integrity.append({'directory':name,'files_checked':len(items),'all_unchanged':True,
                      'archive_sha256':hashlib.sha256(folder.with_suffix('.zip').read_bytes()).hexdigest()})
(H/'prior_round_integrity.json').write_text(json.dumps(integrity,indent=2),encoding='utf-8')
print(json.dumps({'report':str(H/'RESULTS_CN.txt'),'prior_files_checked':sum(r['files_checked'] for r in integrity)}))
