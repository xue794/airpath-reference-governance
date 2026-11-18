"""Package the Section 3 numerical certificate and reader-facing evidence."""
import csv,json,zipfile,hashlib,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from model import *
from audit import run

def tab(cols,rows):
    def row(tag,line):return '<tr>'+''.join(f'<{tag}>{html.escape(str(x))}</{tag}>' for x in line)+'</tr>'
    return '<table><thead>'+row('th',cols)+'</thead><tbody>'+''.join(row('td',r) for r in rows)+'</tbody></table>'
def unit(n):
    if n=='phiMovingBackup':return '椭球二次型/s'
    if n in ['kappa','cosAlpha']:return '1'
    if n in ['Malpha','fbMa']:return 'N/rad'
    if n in ['forceNorm','normal','bankPlus','bankMinus']:return 'm/s²'
    if n=='fbV':return 'm/s'
    if n in ['fbLo','fbHi','rawTmin','rawTmax','payLower','payUpper']:return 'N'
    phi=n.startswith('phi');psi=n.startswith('psi')
    if 'nuV' in n:return 'm/s³' if phi else 'm/s²'
    if 'nuG' in n or 'nuC' in n:return 'rad/s²' if phi else 'rad/s'
    if any(z in n for z in ['Mu','AG']):return 'm/s³' if phi else 'm/s²'
    if 'V' in n:return 'm/s³' if phi else ('m/s²' if psi else 'm/s')
    if n.startswith(('phiG','hG','psiG')):return 'rad/s²' if phi else ('rad/s' if psi else 'rad')
    return 'N/s' if phi else 'N'
def plot():
    f,axes=plt.subplots(1,2,figsize=(11.0,3.45),gridspec_kw={'width_ratios':[1.24,1]})
    zh=FontProperties(fname='C:/Windows/Fonts/msyh.ttc')
    t=np.linspace(0,55,551);v=np.array([center(x,22)[0] for x in t]);h=geometry(.01)['projection'][0]
    ax=axes[0];ax.axhspan(R['Vb']-.606,R['Vb']+.606,color='#e7e8e9')
    ax.fill_between(t,v-h,v+h,color='#96c9bf');ax.plot(t,v,color='#11655e',lw=1.6)
    ax.axvline(35,color='#454e54',ls=':',lw=1);ax.text(36,32.2,'恢复尾段',fontproperties=zh,fontsize=9)
    ax.set(xlim=(0,55),ylim=(30.4,39.2));ax.set_title('移动截面的空速投影',fontproperties=zh)
    ax.set_ylabel('参考空速 (m/s)',fontproperties=zh)
    ax=axes[1];s=np.genfromtxt(OUT/'22_0.01_directional_focused_slabs.csv',delimiter=',',names=True)
    for key,label,col in [('payUpper','上端损失上界','#9b4b1d'),('payLower','下端损失上界','#176778')]:
        ax.stairs(10-s[key],np.r_[s['start'],s['end'][-1]],baseline=None,label=label,color=col,lw=1.7)
    ax.axhline(10,ls='--',color='#363b40',lw=1,label='10 N 预留')
    ax.set(xlim=(14,20),ylim=(6,10.25));ax.set_ylabel('力损失 (N)',fontproperties=zh)
    ax.set_title('关键窗内逐片的保守上界',fontproperties=zh);ax.legend(prop=zh,fontsize=8,frameon=False,loc='lower left')
    for a in axes:
        a.set_xlabel('时间 (s)',fontproperties=zh);a.spines[['top','right']].set_visible(False);a.grid(axis='y',alpha=.14)
    f.tight_layout();f.savefig(OUT/'certificate.png',dpi=155);plt.close(f)
def main():
    run();plot()
    dec=json.loads((OUT/'decision.json').read_bytes());rows=[]
    for mode in dec['whole_tube_bounds']:
        reps=json.loads((OUT/f'22_0.01_{mode}_review.json').read_bytes())
        for k,value in dec['whole_tube_bounds'][mode].items():
            w=min((z['worsts'][k] for z in reps.values()),key=lambda v:v['lower'])
            rows.append([mode,k,unit(k),value,w['slab'],w['path'],*w['time_and_offsets_low'],*w['time_and_offsets_high']])
    # Exact rational rate 3/25 and frozen binary64 rho. Widen displayed lower bound.
    rows.append(['nominal_focused','phiMovingBackup',unit('phiMovingBackup'),float(np.nextafter(.0012,0)),'all','analytic',0,*['']*5,55,*['']*5])
    with (OUT/'all_bounds.csv').open('w',encoding='utf-8-sig',newline='') as fp:
        writer=csv.writer(fp);writer.writerow(['proof','condition','unit','verified_lower_bound','slab','leaf_path',*[p+k for p in ['low_','high_'] for k in ['t','dV','dgamma','dnuV','dnuGamma','dnuChi']]]);writer.writerows(rows)
    show=[['0–55 s','完整五维截面','状态约束、HOCBF 初始化、根单调性、17 条备份输入导数行、第 18 行', '通过'],
          ['14–20 s','分轴误差盒','10 N 分别覆盖上、下端最坏有符号损失，反馈根可用，法向与滚转几何通过','通过'],
          ['14–20 s','分轴误差盒及 θF+θR ≤ 0.70°','投影系数为正，原始反馈推力在 (0,60) N 内','条件成立'],
          ['14–20 s','统一误差盒','计入额外参考裕度时反馈括根和几何通过；10 N 单独支付及原始反馈推力下限各有反例','部分成立']]
    css='body{font:16px/1.7 "Source Sans 3","Microsoft YaHei",sans-serif;color:#173139;max-width:1050px;margin:32px auto;padding:0 25px}h1{font-size:30px;line-height:1.35}h2{font-size:21px;margin:30px 0 6px}p{margin:8px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:7px;text-align:left;border-bottom:1px solid #cad7d7}th{background:#edf4f3}img{max-width:100%;height:auto}a{color:#0c625c}.lead{border-left:4px solid #16736a;padding:7px 14px;background:#eef6f3}small{color:#405b61}'
    report=['<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><title>第 3 节 · 减速机动的可恢复参考管</title><style>'+css+'</style><main>',
      '<h1>22 s 减速中心的连续恢复与力裕度证书</h1><small>Aerospace 大修研究　·　第 3 节　·　2026 年 11 月 8 日</small>',
      '<p class="lead"><b>首选候选通过区域验证，可进入直线开发工况 D0。</b>在名义参考模型上，从原配平到低速保持、再到配平的整个移动管均有可行备份输入。在最关心的 14–20 s，原方案冻结的分轴误差盒允许把反馈括根归因于两端各 10 N 的预留。本阶段没有进行新的六自由度闭环试验。</p>',
      tab(['时间范围','验证域','所核条件','结果'],show),'<img src="results/certificate.png" alt="移动空速截面与逐时间片的反馈力损失上界"><p><small>左图灰带为旧椭球的空速投影，绿带为新管的投影。右图画的是各时间片的区间上界。图中的参考空速与力损失界均不代表实际飞机的飞行记录。</small></p>',
      '<h2>从配平到低速及终端</h2><p>原合法初值的五维误差为零，h(0)=0.01。3–25 s 采用五次平滑中心从 38.1 降至 31.25 m/s，25–35 s 保持；35–55 s 的尾段回到原配平。移动截面的空速投影半宽约 0.191628 m/s。在 25–35 s 整个截面位于约 31.058–31.442 m/s，远离旧椭球 38.1±0.606 m/s 的投影。35 s 是试验计分终点；尾段提供恢复路径。22 s 是备份中心用时，拟测试的外部请求仍用 8 或 10 s。</p>',
      '<p>沿用五维误差 e=(V−Vc, γ, νV−νVc, νγ, νχ)，冻结二进制矩阵 P 与 K，将备份输入取为 uc−Ke。中心满足位置与速率的一致性，故 ė=Aᴮe。对存储浮点数做精确有理运算，可证 −AᴮᵀP−PAᴮ ⪰ 0.12 P。沿备份的第 18 行满足 ḣ+0.4h ≥ 0.12 eᵀPe + 0.4(0.01−eᵀPe) ≥ 0.0012。在线接口在常数项中加入中心的时间导数。各相位的中心值、速度与加速度连续；55 s 中心在原配平且速率和加速度为零，整截面 eᵀPe≤0.01 包含于原 eᵀPe≤0.10，备份输入与原静态策略一致。继续应用该备份时，小截面自身也可保持不变。</p>',
      '<h2>参考递归可行性及其对象</h2><p>在每个时间和截面状态，解析备份给出同时满足运行行的输入，边界内向。采用满足时变第 18 行的连续时间参考更新并从允许初值出发，即可保持在管中，下一时刻仍能使用可行备份。适用前提包括冻结模型、有效上游坐标、名义参考动力学和一致的时间相位。证书使用二进制模型常数及向外区间运算。它不包含采样保持误差、硬件时延或任意请求下实际飞机的跟踪证明。</p>',
      '<p>全域的参考推力比 3 N 下限至少高 2.31×10⁻⁶ N；已扣除 10 N 的上端点参考状态余量下界超过 2.1667 N，备份上端点导数行下界超过 0.0368 N/s。计算在各叶首次严格通过时停止，因此很小的正下界只能用于判号，不能作为可承受的模型误差。逐行表保存最保守叶的时间区间和五个偏移区间。</p>',
      '<h2>两种误差盒给出不同的力预算</h2><p>分轴半宽为空速 0.125 m/s、爬升角 0.45°、最短航向差 0.35°。对整个 [14,20] s × 移动截面 × 误差盒，最坏下端损失不超过 8.0799 N，上端损失不超过 9.9931 N。两端 10 N 预留的已证剩余分别超过 1.9201 N、0.00694 N。结合参考端自身盈余，反馈原始上端力差下界超过 6.2312 N，下端超过 195.71 N。前两数核查预留的作用，后两数核查实际可括根条件。</p>',
      '<p>分轴盒上另证法向及滚转几何。若姿态误差之和 θF+θR≤0.70°，还得到 κ>0 及原始反馈推力 3.18×10⁻⁶ N &lt; Traw &lt; 28.051 N。这组数值属于条件联合证书；真实误差能否留在盒内、姿态条件能否满足、执行器是否触限，需要新闭环逐时核查。</p>',
      '<p>统一盒将空速半宽放大到约 0.29924 m/s，两角半宽均为 0.45°。14 s 的中心上存在允许误差，使上端损失超过 11.9432 N；完整区域计入额外名义裕度后，反馈上端力差仍大于 2.7790 N。在同一时刻另有严格管内点，取允许误差和零姿态偏差，原始反馈推力约为 −2.7687 N。该反例排除统一盒上的 0–60 N 联合保证，不能解读为已发生实际负推力飞行。</p>',
      '<h2>复核与下一步</h2><p>直接区间细分产生未决叶，随后用均值定理、椭球支撑函数及隐式根导数收紧，并把括根和附加推力条件分别验证。包括初次未决搜索在内，名义全管计算累计 987,570 个节点，分轴盒累计 920,562 个节点，统一盒累计 825,934 个节点，分别低于预定 10⁶ 节点及 4 核时限额。旧尝试的局部下界没有用于填补缺失区域。</p>',
      '<p>复核程序从时间片的根开始重建每棵二叉树，对全部最终叶检查覆盖、剪枝或正下界，并用两种三角函数区间求值重算。两次复核共享力学公式。另用 MATLAB 原始约束和新第 18 行核对 1,271 个中心及截面状态；备份输入下 1–17 行与归档最大差约 3.77×10⁻¹³，第 18 行与独立 Lyapunov 恒等式最大差约 1.05×10⁻¹⁷。</p>',
      '<p>下一步先运行 D0（无风、直线、31 m/s、8 s 减速），核对 16 s 前的历史共同前缀、管成员、端点行干预、实际误差盒覆盖与减速效果，再按方案决定是否开放 C1–C3。当前证书给出选择这组工况的依据；相对推力限制或等预留攻角保护的实际闭环增益仍待配对验证。</p>',
      '<p><a href="results/all_bounds.csv">逐行数值、单位与所在叶</a>　<a href="results/decision.json">机器可读判定</a>　<a href="results/raw_thrust_0.01.json">统一盒反例</a>　<a href="README_复算.txt">复算说明</a></p></main></html>']
    (HERE/'第3节_移动参考管区域验证.html').write_text('\n'.join(report),encoding='utf8')
    files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='release_manifest.json']
    manifest={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    fp=HERE/'release_manifest.json';save(fp,manifest)
    zpath=ROOT/'Aerospace_有限时长减速机动_第3节区域验证_v1_0.zip'
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in [*files,fp]:z.write(p,p.relative_to(HERE).as_posix())
    with zipfile.ZipFile(zpath) as z:
        assert z.testzip() is None and all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in manifest.items())
    print(zpath.stat().st_size)
if __name__=='__main__':main()
