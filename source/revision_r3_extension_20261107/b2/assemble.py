"""Build the checked case tables, figure and Chinese B2 research report."""
from collections import Counter
from fractions import Fraction as F
import csv,gzip,json,math
from html import escape
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from coverage import H,ROOT,A3,writetable
from search import writej

def csvrows(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def tab(labels,rows):
    make=lambda tag,vs:"".join(f"<{tag}>{escape(str(v))}</{tag}>" for v in vs)
    return '<div class="scroll"><table><thead><tr>'+make("th",labels)+'</tr></thead><tbody>'+''.join(
        "<tr>"+make("td",r)+"</tr>" for r in rows)+'</tbody></table></div>'
def figure(rows):
    fig,axs=plt.subplots(1,2,figsize=(10.4,4.1),width_ratios=[1.15,1],layout="constrained")
    colors=["#176781","#aa501c"]
    marker={"pass":"o","witness":"x","unresolved":"^"}
    for panel,ax in enumerate(axs):
        subset=[r for r in rows if (r["stage"]=="grid" if panel==0 else .18<=r["eps"]<=.222)]
        ax.axhline(10,color="#333",ls="--",lw=1)
        for j,key in enumerate(["U_lower_N","U_upper_N"]):
            for r in subset:
                wkey=["witness_lower_N","witness_upper_N"][j]
                state="pass" if r[key]<10 else ("witness" if r[wkey]!="" else "unresolved")
                ax.scatter(r["eps"],r[key],s=30,marker=marker[state],color=colors[j],
                    facecolors="none" if state=="unresolved" else colors[j],zorder=3)
        ax.set_xlabel("Joint error half-width ε (m/s)")
        ax.grid(axis="y",alpha=.18);ax.spines[["right","top"]].set_visible(False)
        ax.set_title(["Prespecified grid","Refinement around the 10 N threshold"][panel],fontsize=11)
    axs[0].set_ylabel("Upper bound on endpoint loss (N)")
    axs[0].set_ylim(-1.5,57);axs[1].set_ylim(6.8,11.7)
    axs[1].axvline(51/256,color="#666",lw=.7,ls=":")
    handles=[Line2D([],[],color=c,marker="s",ls="",label=k) for c,k in zip(colors,["Lower endpoint","Upper endpoint"])]
    handles.extend(Line2D([],[],color="#333",marker=m,ls="",label=l,markerfacecolor="none" if m=="^" else "#333")
                   for m,l in [("o","Bound < 10 N"),("^","Unresolved"),("x","Witness > 10 N")])
    fig.legend(handles=handles,loc="outside lower center",ncol=5,frameon=False,fontsize=8)
    fig.savefig(H/"reserve_curve.png",dpi=200);plt.rcParams["svg.hashsalt"]="B2-fixed-eps"
    fig.savefig(H/"reserve_curve.svg",metadata={"Date":None});plt.close(fig)

def assemble():
    j=json.loads((H/"search.json").read_bytes());qa=json.loads((H/"tree_review.json").read_bytes())
    wi={r["epsilon"]:r["witnesses"] for r in json.loads((H/"counterexamples.json").read_bytes())}
    rows=[]
    for i,r in enumerate(j["rows"]):
        w=wi.get(r["epsilon"],[]);status="pass" if r["passed"] else ("witness" if w else "unresolved")
        rows.append(dict(epsilon=r["epsilon"],eps=r["eps"],stage="repeat" if r["repeat"] else ("grid" if i<8 else "bisect"),
            status=status,U_lower_N=r["max_lower_pN"]/1e12,U_upper_N=r["max_upper_pN"]/1e12,
            witness_lower_N=max([v["loss_lower_N"] for v in w if v["side"]=="lower"],default=""),
            witness_upper_N=max([v["loss_lower_N"] for v in w if v["side"]=="upper"],default=""),
            visited=r["visited"],pruned=r["pruned"],accepted=r["accepted"],unresolved_leaves=r["unresolved"],tree=r["tree"]))
    writetable(H/"radius_results.csv",rows);figure(rows[:-1])
    run=csvrows(H/"coverage_by_run.csv");ev=csvrows(H/"event_coverage.csv");groups={}
    for name,test in [("保持",lambda r:r["group"]=="WP4" and r["trajectory"].endswith("hold")),
        ("时变风",lambda r:r["group"]=="WP4" and r["trajectory"].endswith("varying")),
        ("含噪声",lambda r:r["group"]=="WP4" and r["trajectory"].endswith("noisy")),
        ("A1 原 QP",lambda r:r["group"]=="A1"),("31 m/s 压力",lambda r:r["group"]=="stress" and r["trajectory"].startswith("V31")),
        ("32 m/s 压力",lambda r:r["group"]=="stress" and r["trajectory"].startswith("V32")),
        ("A3 配置记录",lambda r:r["group"]=="A3")]:
        rr=[r for r in run if test(r)]
        groups[name]={k:sum(int(r[k]) for r in rr) for k in ["n","inside","boundary","outside","B1_016","covered","G","J_air_box","J","KB_extra_1"]}
        groups[name]["records"]=len(rr)
    writetable(H/"coverage_groups.csv",[dict(group=n,**v) for n,v in groups.items()])
    evt={}
    for name,ee in [("initial",[r for r in ev if float(r["time_s"])==0]),("inflight",[r for r in ev if float(r["time_s"])>0])]:
        evt[name]=dict(n=len(ee),**{k:sum(int(r[k]) for r in ee) for k in ["covered","G","J","B1_016","KB_extra_1","old_local_J"]})
    critical=[]
    for name in ["V31_N10","V32_N10"]:
        with gzip.open(H/"coverage"/("stress_"+name+".csv.gz"),"rt") as f:record=list(csv.DictReader(f))
        smallest=min(record,key=lambda r:float(r["fb_upper_N"]))
        critical.append(dict(case=name,end_covered_s=max(float(r["t"]) for r in record if int(r["E"])),
            minimum_upper_s=float(smallest["t"]),minimum_upper_N=float(smallest["fb_upper_N"]),at_min_covered=int(smallest["E"])))
    done=rows[-1];last=F(done["epsilon"])
    uncertain=min(F(r["epsilon"]) for r in rows if r["status"]=="unresolved")
    failed=min(F(r["epsilon"]) for r in rows if r["status"]=="witness")
    assert last<uncertain<failed and done["eps"]==51/256 and len(ev)==21
    summary=dict(radius=str(last),epsilon_mps=float(last),angle_deg=float(last)*180/(38.1*math.pi),
       next_unresolved=str(uncertain),first_strict_witness=str(failed),lower_loss_bound_N=done["U_lower_N"],
       upper_loss_bound_N=done["U_upper_N"],verified_distinct_leaves=qa["checked_distinct_passing_tree_leaves"],
       groups=groups,event_groups=evt,pressure_minima=critical)
    writej(H/"b2_verified.json",summary)
    h='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>B2 · 原备份椭球上的端点预留与轨迹覆盖</title><style>
body{max-width:1080px;margin:auto;padding:34px min(5vw,54px) 66px;font:16px/1.72 "Microsoft YaHei",sans-serif;color:#192e38}
h1{font-size:29px;line-height:1.3}h2{font-size:21px;color:#174e62;margin-top:28px}p{margin:8px 0}
table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}td,th{text-align:right;white-space:nowrap;border-bottom:1px solid #d2dfe0;padding:5px 9px}
th{background:#eaf1f1}td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}.eq{background:#eef3f4;padding:12px}img{width:100%;height:auto}a{color:#175a71}small{color:#42535b}
@media print{body{padding:0}.scroll{overflow:visible}tr{break-inside:avoid}}</style><main><h1>B2　10 N 预留的区域界与闭环覆盖</h1>
<p><strong>原 ρ<sub>B</sub> = 0.10 椭球上的充分联合半宽从 0.160 提高到 0.19921875 m/s，约扩大 24.5%。</strong>只使用参考上下端点各 10 N 的预留，就能在该误差盒内保证名义反馈反演在原 0°–9° 内有唯一根。盒内允许空速误差达到 ±0.19921875 m/s，同时允许爬升角和最短航向误差各达到 ±0.19921875/38.1 rad（约 ±0.2996°）。结论以当前参考属于椭球、误差属于盒为条件。</p>
<p>在 A2 的 21 个确认作用点中，5 个初始化点落在该区域，9 个飞行中作用点全部在误差盒外。A3 样本 22 的 0.182733 N 闭环增益属于直接仿真证据。该例飞行中首次作用时的联合误差为 0.337374 m/s（向上取整），超过本次证得的半宽。</p>
<h2>搜索、反例和数值界限</h2>
<p>最终通过的树访问 91,803 个节点，接受 27,950 个终端盒并剪去 17,952 个椭球外盒。下、上端点的力损失上界分别为 7.546821572388 N、9.999966660617 N；10 N 预留下的剩余间隙至少为 2.453178427612 N、0.000033339383 N。上端点余量很小。紧邻的 409/2048 = 0.19970703125 m/s 尚有 8,570 个未决终端盒。7/32 = 0.21875 m/s 存在域内上端点损失超过 10.76985 N 的严格见证。由此，10 N 足以支付全域损失的临界半宽位于已证的 0.19921875 与严格反例的 0.21875 m/s 之间，当前计算没有确定其精确值。</p>
<img src="reserve_curve.png" alt="两侧反馈力损失的区域上界及10 N门槛。左图为冻结网格，右图为细化结果。"><p><small>蓝色为下端，橙色为上端。实心点表示该侧区域上界小于 10 N，空心三角为该侧预算内未决，叉号表示该侧存在独立验证的严格反例。纵轴显示区域上界；用于反证的一点实值另列。不同半宽采用不同划分，上界之间不作插值。</small></p>'''
    def upfmt(x):return f"{math.ceil(x*1e6)/1e6:.6f}" if x else "0"
    labels={"pass":"已证","unresolved":"未决","witness":"严格反例"}
    h+=tab(["半宽 (m/s)","阶段","下端上界 (N)","上端上界 (N)","状态","上端见证损失下界 (N)","未决盒"],
        [[str(r["eps"]),{"grid":"原网格","bisect":"细化"}[r["stage"]],upfmt(r["U_lower_N"]),upfmt(r["U_upper_N"]),
          labels[r["status"]],f'{math.floor(r["witness_upper_N"]*1e5)/1e5:.5f}' if r["witness_upper_N"]!="" else "—",r["unresolved_leaves"]] for r in rows[:-1]])
    h+='''<h2>原始轨迹覆盖</h2>
<p>每一行从该组所有保存点开始计数。先核对名义力核和矩阵，再要求参考可判定属于原椭球、联合误差落在获证盒内，得到 E。G 进一步核查法向请求与坡度角条件。J 还要求姿态误差和 ≤0.70°、正推力投影、整盒的原始反馈推力严格落在 (0,60) N，并核实当拍的原始请求。实际执行器推力和参考推力另列在逐轨迹台账中。边界栏指按相邻浮点数外包后无法判定的椭球成员，已从证书覆盖数排除。</p>'''
    h+=tab(["样本组","轨迹数","保存点","椭球内","边界未定","椭球外","B1 0.160","B2 E","G","J","额外间隙 1 m/s 的 E"],
        [[n,*[f'{v[k]:,}' for k in ["records","n","inside","boundary","outside","B1_016","covered","G","J","KB_extra_1"]]] for n,v in groups.items()])
    h+='''<p>保持、时变风和含噪声各含 Scalar、Matched-QP、Transparent。A1 是按原上端点支持筛选的 17 条完成轨迹。其 4,385 个可判在椭球外的保存点均在原 A1 数值容差 1+10⁻⁸ 内，严格的区域成员与积分验收容差分别统计。A3 的 18 条是三例、三配置、两档最大步长的记录，其中三条 N10 原步长来自 A1；这些行不可合并为独立试验成功率。</p>
<p>压力组采用已存的细步长记录。31 m/s 包含 T、N0、N10、A0、Am，32 m/s 包含 T、N0、N10、Am；总数包含提前终止时实际保留的末点。两条 N10 轨迹各有 701 点，区域证书分别只覆盖前 112、114 点，最后覆盖时刻为 5.55、5.65 s。其最小反馈上端点分别出现在 18.70、18.40 s，均在证书覆盖时间之外。逐点直接检查与先前局部误差盒结果保留各自口径。</p>
<h2>端点行实际作用的 21 点</h2>'''
    h+=tab(["时段","A2 确认点","B1 0.160","B2 括根","B2 联合推力","额外间隙 1 m/s 括根","原逐点局部 J"],
       [[name,*[v[k] for k in ["n","B1_016","covered","J","KB_extra_1","old_local_J"]]] for name,v in zip(["初始化","飞行中"],evt.values())])
    h+='''<p>本次覆盖的五个作用点为 80、91、164、386、483 的初始化点。它们均通过 G 且当拍原始推力处于范围内，但姿态误差和均超过 0.70°，故此处区域 J 的交集为零。A1 的局部 J 在另外三个初始化点通过，使用各自参考状态处的局部误差盒。若允许使用原备份集额外的参考力间隙，B1 的 1 m/s 括根结果覆盖 16/21 点，其中 4/9 为飞行中点。这个较宽保证依赖整个备份构造。</p>
<h2>与 A3 的闭环代价并列</h2><p>下表复用原 120 s 或原终止规则的 A3 细步结果。反馈指标在首次 A2 确认时刻至共同终点取最小值，正差表示 N10 保留的反馈上端点间隙更大。RMS 使用自 0 s 起的全共同网格。三种配置在这三例均完成 120 s，推力请求裁剪计数为零。原 500 对试验仍为两法各 431/500 完成；普通工况六条轨迹的 12,006 个保存点中，原同状态端点删行指令零差异。</p>'''
    con=csvrows(A3/"contrasts.csv");sr=csvrows(A3/"run_summary.csv")
    h+=tab(["编号","窗口起点 (s)","N10−D16 (N)","N10−Am (N)","全程 RMS N10 / D16 / Am (m)"],
        [[i,next(r["start_s"] for r in con if r["id"]==i),*[f'{float(next(r["fine_delta_N"] for r in con if r["contrast"]==f"N10-{a}" and r["id"]==i)):.7g}' for a in ["D16","Am"]],
          " / ".join(f'{float(next(r["full_common_RMS_m"] for r in sr if r["id"]==i and r["arm"]==a and r["step"]=="fine")):.6f}' for a in ["N10","D16","Am"])] for i in ["22","32","33"]])
    h+='''<h2>区间方法和复验</h2>
<p>沿用 B1 已复核的 MATLAB 二进制名义参数和 ρ<sub>B</sub>=0.10 投影外盒。B1 的运行参数有理证明及旧／新参数共同外包已重新核对参考间隙。对 ε≤1，B1 给出的 V&gt;36.4940 m/s、|γ|&lt;0.05357 rad、全角域 ∂M/∂a&gt;954.55 N/rad 覆盖所有子盒。旧十进制证书的额外裕度没有用于 10 N 的通过判定。</p>
<p>为在区域盒中保留相消，记 p=(A,H) 为竖直与转弯方向的单位质量法向请求。对参考和反馈两端分别计算 p，使用 ΔN=m w·(p<sub>c</sub>−p<sub>f</sub>)，其中 w=(p<sub>c</sub>+p<sub>f</sub>)/(‖p<sub>c</sub>‖+‖p<sub>f</sub>‖)。w 是两端单位向量的凸组合。用 A&gt;0 及 H/A 的包络求其两个分量界，再将 ΔN 与端点 ΔM 中同一误差坐标的系数组合。正弦、余弦差用参考与反馈角之间的导数包络。两侧损失分别为 ΔM(0°)−ΔN 和 ΔN−ΔM(9°)。</p>
<p>主实现采用 56 位小数固定点外向运算，超越函数单点用相邻有理 Taylor 部分和取界。椭球二次型配方后，只有整盒下界严格超出 ρ<sub>B</sub> 才剪枝。按归一化宽度平局顺序循环选维，广度优先、低半盒在前。每次生成两个子盒，100,000 预算下最大完整二叉树有 99,999 个节点；所有待细分终端都保存了其上界。区域总上界取所有非剪枝终端（含未决）的最大值。全域已有正半宽通过，未触发 0.025 子域。独立的 65 位十进制外向实现复验所有通过树的 52,366 个终端盒，校验全体树的拓扑和最坏盒中心，最终通过树重跑逐字节相同。严格反例先用数值方法定位，再以有理坐标用两套外向实现核对域、误差盒和损失符号。</p>
<p><b>论文可写入的结论。</b>可在原认证域内由反馈误差计算参考预留，并用严格见证界定 10 N 规则的适用范围。普通保持工况获得全保存网格覆盖。随机样本的飞行中作用点、压力段的最紧时刻主要由直接闭环证据支持。误差持续处于盒内的闭环保证、全任务安全和任务意义上的跟踪改善仍需各自证据。</p>
<p><a href="radius_results.csv">搜索结果</a> · <a href="coverage_by_run.csv">53 条逐轨迹结果</a> · <a href="event_coverage.csv">21 点交集</a> · <a href="counterexamples.json">反例坐标</a> · <a href="tree_review.json">独立复核</a> · <a href="README_复算说明.txt">复算说明</a></p></main></html>'''
    (H/"B2_区域预留与轨迹覆盖.html").write_text(h,encoding="utf-8")
if __name__=="__main__":assemble()
