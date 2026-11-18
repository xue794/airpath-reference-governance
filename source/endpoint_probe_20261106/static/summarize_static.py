from pathlib import Path
import csv, json

root = Path(__file__).resolve().parent
summary = json.loads((root/'static_summary.json').read_text())
interval = json.loads((root/'interval_scope.json').read_text())
rows = list(csv.DictReader((root/'static_comparisons.csv').open(encoding='utf-8-sig')))
selected = [r for r in rows if float(r['nuV']) == 0 and float(r['nuGamma_deg']) == 0
            and float(r['nuChi_deg']) == 0 and float(r['V']) in [31, 38.1]
            and r['target'] in ['pullup', 'turn30']]
(root/'selected_witnesses.json').write_text(json.dumps(selected, ensure_ascii=False, indent=2), encoding='utf-8')
text = '''第二轮静态端点约束探测，2026-11-06

问题与比较边界
仅推力基线已经使用与完整方案相同的耦合隐式反演根与参考推力Tff。
因此本实验隔离的是端点可行性导数约束的附加作用，不检验耦合反演本身
相对于独立升力/推力近似的优越性，也不与全部既有包线方法作比较。

三组约束为：
A：1:13和16:17（运动学/几何+同耦合Tff的3–57N推力约束）；
B：1:17（A再加两条端点约束，DeltaN=10N）；
C：1:17但DeltaN=0（区分根括号边界与人为10N预留）。
三组使用相同历史QP、权重、原始请求和同一参考状态；均不包含KB。
由于不含KB，原论文的递归可行性保证不适用于这些比较。

为避免初始条件不公平，只有同时满足以下条件的状态进入比较：
常规速度、爬升、率、正法向及银行角状态条件；四个HOCBF psi1非负；
当前0–9deg反演根存在；Tff在3–57N；两端点10N预留均满足。
有限格：V=[30,31,32,34,36,38.1]m/s，gamma=chi=0，
nuV=[-0.5,0,0.5]m/s²，nuGamma=[0,0.5,1,2,3]deg/s，
nuChi=[0,2,4,6,8]deg/s。7个请求包含爬升、下降、转向、低速爬升等。
该网格是公开的探索性扫描，其命中比例不是实际运行概率。

结果
450个状态中277个通过全部上述初始化检查，得到1939次同状态请求比较。
三组QP在全部1939次比较中均得到可验证解，最大归一化KKT残差6.62e-16。
仅推力解有629次违反10N端点预留的导数条件；其中583次在DeltaN=0时
仍违反原始端点导数条件。故这种作用不能全部归因于10N预留。
但导数行违反不是状态违反，更不是未来必然失根；动态任务另行验证。

简单代表例
当前x=(V=31m/s,gamma=0,chi=0,nuV=nuGamma=nuChi=0)，
alpha=8.094339deg，Tff=24.782992N，原始上端点裕度19.294158N，
扣除10N预留后剩余9.294158N。所有比较组的初始条件均满足。
请求c=(31m/s,8deg,0)：
- 仅推力解u=(-0.0181892,0.2004054,0)，原始上端点导数条件=-120.542149N/s；
- 完整10N预留解u=(0.1957508,0.0143262,0)，上端点导数条件约为0；
- 完整0N预留解u=(0.1808006,0.0285879,0)，原始上端点导数条件约为0。
u的单位依次为m/s³、rad/s²、rad/s²。
这一例证明在当前根存在且推力约束可行时，端点约束仍能改变可接受更新。
它本身没有证明这一个31m/s状态能由原始KB约束方案达到。

反演区间与物理诊断区间
代码的实际控制反演区间是0–9deg；alphaPhysicalMin/Max=-2/10deg
是研究中采用的另一个诊断范围，并非由本实验认证的真实硬件包线。
V31m/s、nuChi=8deg/s、其他率为0的状态在0–9deg内无根，
但在仅作反事实诊断的[-2,10]deg范围内有根：alpha=9.060362deg，
Tff=30.059330N。故0–9deg求根失败不能直接称为物理攻角越界。
V31m/s、nuGamma=4deg/s例在扩至10deg后仍不括根，也仍然只是
所给标称反演、所选诊断区间下的不可实现性，不能泛化到所有控制架构。

文件与复现
run_static_endpoint_probe.m：完整扫描；MATLAB R2018a，自动引用第一轮冻结核心。
static_states.csv：所有450个状态，含被筛除点及原因相关量。
static_comparisons.csv：全部1939次比较，三组控制向量、裕度与KKT残差。
static_endpoint_probe.mat：完整参数、数据与汇总。
audit_interval_scope.m / interval_scope.json：反演区间反事实诊断。
selected_witnesses.json：简单具名展示实例，均从完整扫描表中挑选。
summarize_static.py：生成本文档和具名实例清单。

建议对论文的有限结论
“仅保证耦合参考推力处于区间内，不足以保证选定反演区间的导数
可行条件；端点行可在推力约束已经存在时独立限制参考更新。”
是否能在合法初态出发的闭环中避免失根、其跟踪代价及是否优于直接
alpha约束，需要动态试验和明确基线继续验证，不能从本静态扫描推断。
'''
(root/'STATIC_RESULTS.txt').write_text(text, encoding='utf-8')
print(json.dumps(summary, indent=2))
