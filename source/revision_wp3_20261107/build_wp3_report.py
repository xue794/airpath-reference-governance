from pathlib import Path
import json,math
from copy import deepcopy
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
H=Path(__file__).resolve().parent
S=json.loads((H/'results/summary_table.json').read_text(encoding='utf-8'))
T=json.loads((H/'results/transfer_summary.json').read_text(encoding='utf-8'))
M=json.loads((H/'results/margin_summary.json').read_text(encoding='utf-8'))
name={'Proposed':'Scalar','MatchedQP':'Matched-QP','Transparent':'Transparent'}
case={'hold':'H 静风保持','varying':'V 时变风','noisy':'N 测量误差'}
def flo(v):return f'{math.floor(v*1e6)/1e6:.6f}'
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
P=[]
def page(title,*blocks):P.append((title,blocks))
def p(x):return ('p',x)
def sub(x):return ('sub',x)
def tbl(h,r,w):return ('table',h,r,w)
def eq(x):return ('eq',x)
def fig(x,width):return ('figure',x,width)

english='''Quantitative local transfer (proposed insertion after Proposition 4)

We quantify the conditional instantaneous transfer of the reference margins on the retained original local examples. For a fixed speed scale of 38.1 m/s, the air-path uncertainty is a componentwise box with speed half-width epsilon and angular half-widths epsilon/38.1. The interval test retains the order of Proposition 4: root-free normal, bank, and endpoint conditions are checked before establishing the unique implicit root and bounding raw feedback thrust through the attitude-dependent projection factor. The original 0.70-degree attitude-sum bound is used as a condition, not as an invariant error guarantee.

All 14,407 saved samples from the seven available original trajectories admit positive joint sufficient half-widths. The geometry-only calculation reaches the imposed 1 m/s search ceiling everywhere and covers all recorded air-path errors; this ceiling is not a maximum admissible radius. Joint inclusion, additionally requiring raw feedback thrust in (0,60) N and the original attitude condition, is 1,201/1,201 in each hold run, 1,690/2,401 and 1,717/2,401 for Scalar and Matched-QP in varying wind, and 1,688/2,401 and 1,712/2,401 with the original deterministic measurement error. Transparent varying gives 1,582/2,401. The finite-sample minimum joint half-widths are approximately 0.240108 and 0.229284 m/s for Scalar and Matched-QP in varying wind, respectively; the exact sufficient values and their units are retained in the data record.

The attitude-sum condition is observed at every saved point. For Scalar and Matched-QP, the remaining varying/noisy samples fail the conservative lower-thrust enclosure over the whole observed-error box. None of the seven directly reconstructed feedback margins is negative at any retained sample. The sufficient-condition gap therefore does not identify a physical violation. These results quantify a local command-transfer condition and its limitations, without establishing whole-ellipsoid uniformity, between-sample guarantees, or invariance of the tracking-error box. They also do not establish a benefit unique to the endpoint representation relative to a matched angle-of-attack constraint.

Response to Reviewer 2 on the unspecified neighborhood (draft)

We have added a quantitative, ordered interval construction to Supplementary S4 and a table for the original local examples. The revision specifies the scaled error coordinates, the conditional attitude bound, all seven checked command margins, the search ceiling, and the positive-certificate rule. Endpoint signs and monotonicity are verified before the root enclosure and thrust calculation. We report both the sufficient half-widths and the recorded-condition inclusion counts, including the uncovered samples. The result is a sampled-reference instantiation of the local command-transfer statement, not a uniform certificate over the full backup set or a closed-loop error-invariance proof. The manuscript's guarantee table has been narrowed accordingly when these additions are integrated.
'''
(H/'manuscript_insert_EN.txt').write_text(english,encoding='utf-8')

page('Aerospace 原局部案例反馈裕度传递核查',
 p('WP3 完成报告  2026年11月7日'),
 p('原局部案例的反馈裕度传递已得到定量结果。7条既有轨迹的14,407个保存点均有正的联合充分半宽；端点、法向加速度和倾侧角条件的误差盒检查覆盖全部保存点。加入原始反馈推力条件后，时变风及测量误差案例的联合条件覆盖约70%，没有形成全段覆盖。'),
 p('未覆盖不等于违约。按原模型直接重建的七项反馈命令裕度在全部保存点都为正；这些结果支持“条件性的局部传递”，不能转成连续时间整机安全、误差管不变性或专用端点约束的独立性能优势。'),
 tbl(['核查项目','结果'],[
 ['复用范围','hold 两条；varying 三条；noisy 两条。新增轨迹积分为0。'],
 ['参数与方法','原十八行、rho=0.10、原增益与风/参考波形，DeltaN=10 N、DeltaT=3 N。'],
 ['关键数值','Scalar varying 1690/2401；Scalar noisy 1688/2401；hold 两方法均1201/1201。'],
 ['核查强度','MATLAB原函数重建；40位区间计算；70位独立角点核查；输入与源码哈希保持一致。'],
 ['论文作用','补齐命题4与S4的具体数值范围；保留未覆盖部分及局部范围限制。'],
 ],[1.15,5.55]),
 p('本阶段保持原研究边界。没有开发恢复备份、改变控制器或扩大认证集合；没有重跑Monte Carlo。暂无独立实际任务指标，报告不设置任务合格线。正文和补充材料原文件均未覆盖，修订文本单独交付。'))

page('一 核查对象与固定条件',
 tbl(['对象','范围或定义'],[
 ['轨迹','H：0–60 s；V/N：0–120 s；统一0.05 s保存网格。两方法H/V/N加Transparent V。'],
 ['身份','7个输入MAT文件SHA256与第一阶段一致；315个代码/支持文件复制并校验。'],
 ['初态与配置','36维保存初态与原初始化函数一致；C与最终名义配置一致。Transparent仅将审计rho从0.20统一为0.10，其选择不使用KB。'],
 ['反馈坐标','从地速扣除真实风，重建空速、航迹角和包裹航向误差。noisy保留原确定性测量误差来源。'],
 ['命令与执行器','Tstar为理想耦合推力；Traw为执行器限幅前原始反馈推力指令；Tact为执行器状态。'],
 ['未复用格','Transparent hold/noisy尚无可识别原记录；不计入本次结论，留给WP4按既定协议补缺。'],
 ],[1.15,5.55]),
 sub('误差盒与姿态条件'),
 eq('box'),
 p('半宽ε的单位为m/s；两个角度误差以弧度计。它是逐分量盒半宽，不是S17中的欧氏半径。缩放后的欧氏球半径ε包含于该盒中，反向包含不能成立。'),
 eq('attitude'),
 p('主计数使用原S4的滤波误差与姿态跟踪误差之和。直接命令姿态距离另列，两者不能暗中替换。0.70°来自原配置，只作为此次数值条件；其未来保持性仍未证明。'))

page('二 先验证无需求根的条件',
 p('反演适用域是原0–9°命令根区间。原始端点裕度用于判断括根；10 N是参考设计预留，不能在反馈原始裕度中再次扣除。固定气动多项式后，端点裕度损失有下列精确差分关系。'),
 eq('margins'),
 p('其中K为固定名义密度与机翼面积乘积的一半，ΔX与ΔN分别是反馈和参考轴向需求、法向合力需求之差。相同参考点与同一误差盒用于三层核查，不将不同方法闭环上的状态混为同一点。'),
 tbl(['层次','需要对整个误差盒成立的条件'],[
 ['E 端点层','两原始端点裕度严格正；原根区间上M的攻角导数下界严格正；空速和航迹坐标正则。'],
 ['G 几何层','E层加上反馈法向加速度大于4 m/s²、两项35°倾侧角仿射裕度严格正。'],
 ['J 联合层','G层加上姿态条件、正投影因子，以及原始反馈推力严格处于0–60 N。'],
 ],[1.1,5.6]),
 p('在所有保存参考点，E与G层都通过预先设置的1 m/s搜索上限，且均包含实际误差。对应两个角度半宽约1.504°。这说明本批轨迹附近的这些几何条件有较充分的局部余量，不能把1 m/s写成真正最大半径。'),
 p('本次量化的是命题4所列的命令裕度，未把实际飞机攻角、所有舵面限制、位置跟踪稳定性或整个KB的统一反馈邻域并入同一证书。'))

page('三 隐式根与推力的有序传递',
 p('先由两端点符号建立根存在，再由整个原反演区间的正导数下界建立唯一性。只有完成这两步，才用参考根近似作为中心，通过残差包围反馈根；没有预先假定要证明存在的根。'),
 eq('root'),
 p('残差区间包含近似中心自身的求根误差；根包络与原0–9°区间取交集。命令姿态使用实际航迹角，因而姿态距离可给出机体轴推力投影因子的几何界。'),
 eq('thrust'),
 p('角度和姿态上界之和还必须小于90°，使投影下界为正。这里限制的是限幅前原始指令，不能把它当成实际执行器推力。0–60 N是原模型推力盒；3–57 N是参考预留，反馈阶段不要求再保持该预留。'),
 sub('数值充分界的含义'),
 p('固定源文件的binary64参数，在40位区间算术内形成组合系数并作外向包围。半宽搜索上限为1 m/s，二分20次，每个正返回值再检验。二分值是通过充分判据的半宽，不是最大可行半径。'),
 p('这是一种对正裕度的数值充分检查，未计算S17/S18中所有Lipschitz常数，也不声称保持原证明中“一半裕度”的具体下界。有限保存点半宽的最小值只对应这些参考点，不能代替整个椭球或采样间轨迹上的统一下界。'))

page('四 联合传递半宽和记录覆盖',
 p('下表最小半宽向下保留6位小数，可作为该组保存参考点共同的充分半宽；每个时间点的半宽通常更大。覆盖按实际误差是否进入该时间点的盒、且姿态和不超过0.70°同时判断。'),
 tbl(['案例','方法','J最小半宽\nm/s','联合覆盖\n保存点','记录比例'],[
 [case[x['case']],name[x['method']],flo(x['J_min_mps']),f"{x['joint_S4_count']}/{x['samples']}",f"{x['joint_S4_percent']:.2f}%"] for x in S
 ],[1.1,1.15,1.35,1.65,1.45]),
 p('所有7条轨迹、14,407个保存点的J层半宽均为正；E/G层在全部点达到搜索上限并包含实际误差。所有点的姿态和也满足0.70°：H最大约0.007716°，V约0.218395°，N约0.212234°。改用直接姿态距离，主覆盖计数完全相同。'),
 tbl(['案例与方法','J半宽中位数\nm/s','J半宽最大值\nm/s','实际误差最大值\nm/s'],[
 [case[x['case']]+' '+name[x['method']],f"{x['J_median_mps']:.6f}",f"{x['J_max_mps']:.6f}",f"{x['max_error_box_mps']:.6f}"] for x in S if x['case']!='hold'
 ],[2.5,1.4,1.4,1.4]),
 p('中位数和最大值为描述性四舍五入值，精确通过值保存在逐点CSV。记录比例不是安全概率、连续时间覆盖率或任务可靠性。不同方法在各自轨迹上计算，不能据此作同状态约束集大小比较。'))

page('五 未覆盖点与实际命令裕度',
 p('Scalar和Matched-QP在V/N中共有2797个未覆盖保存点；它们全部是对整个实际误差盒的推力下界未能证实为正。端点、法向、倾侧角、投影因子和0.70°姿态条件均不是这些点的缺口。Transparent V的819个未覆盖点中，808个对应推力下界，11个对应上界。'),
 p('这是当前充分包围的缺口。它允许比实际误差更不利的符号组合，也包含区间相关性损失；仅凭失败下界不能证明盒内必有违约点。无需据此改控制器或扩大本篇研究。'),
 tbl(['案例与方法','参考推力最小\nN','原始反馈推力\n最小至最大 N','位置RMS\nm'],[
 [case[x['case']]+' '+name[x['method']],f"{x['observed_min_reference_T_N']:.4f}",f"{x['observed_min_Traw_N']:.4f}–{x['observed_max_Traw_N']:.4f}",f"{x['RMS_grid_m']:.4f}"] for x in S
 ],[2.5,1.2,1.8,1.2]),
 p('直接重建的七项反馈命令裕度全部为正，原始反馈推力在这批记录中甚至均未消耗至3 N下预留。建模执行器及已检查研究状态界限也未发现保存点越界。上述均是记录事实，不是连续时间整机保证。'),
 p('Transparent V仍有155/2401点参考推力低于3 N，2339/2401点在KB外，但其原始反馈推力最小约11.20 N且跟踪更准。因此，较大的充分半宽可以作为局部可核查性证据，不能自动解释为这组跟踪代价换来了必要物理保护。'))

page('六 原时变风 Scalar 的传递范围',
 fig('varying_scalar_transfer.png',6.65),
 p('图1  上图蓝线为J层充分半宽，黑线为实际缩放气路误差；浅红色仅标记联合条件未覆盖的保存点附近网格区段。灰虚线是G层通过的搜索上限，不能解释为最大半径。下图为原S4姿态误差和及直接姿态距离，两者均小于固定条件0.70°。'),
 p('曲线连接只帮助阅读保存值，不能据此断言区间内每一时刻通过或失败。完整三方法图、推力分层图和所有逐点数值另外提供；这里选用原V案例作代表，不按结果更换波形。'),
 p('未覆盖区段的首末保存时刻与点数保存在uncovered_sample_segments.json。它们不是经过事件定位的连续时间失效区间，也不作为积分停止或飞机失稳的判据。'))

page('七 原命题与修订文字应如何衔接',
 tbl(['问题','本次可以写入的结果','仍不能写的内容'],[
 ['命题4邻域未量化','7种case/method记录上的三层半宽、全部条件和覆盖计数','整个KB的统一反馈半径；全时域误差管不变性'],
 ['0.70°阈值依据','以原阈值为条件完成联合区间检查，并报告实际姿态和','该阈值已成为控制误差的先验不变上界'],
 ['反演根逻辑','端点符号先成立，再导数和残差包根，再传递推力','预先假定反馈根存在后再证明其存在'],
 ['未覆盖保存点','公开推力充分下界未通过及正的实际重建裕度','把未认证当成违例；把采样无违例当成连续证明'],
 ['审稿人3的贡献','具体耦合接口具有可计算的条件传递实例','端点优于强匹配攻角；专用约束独有的工程性能收益'],
 ],[1.1,2.8,2.8]),
 sub('建议整合位置'),
 p('正文命题4后增加误差坐标、三层计算定义与一张精简结果表；S4增加精确差分、区间判据及保证范围说明。结果节沿用V案例呈现覆盖缺口，讨论段将“有可计算局部余量”与“实际安全收益”分开。'),
 p('本次结果不改变主线为新恢复备份研究。对审稿人2，它补上了数值邻域缺口的一部分；对审稿人3，它增加可复核技术细节，但也继续表明原局部案例中的端点并非限制传递范围的主要环节。不能据此认定核心贡献已得到充分工程验证。'),
 p('WP3按既定有限范围完成并停止。下一项仍是WP4按冻结协议补Transparent hold/noisy记录，并判断是否需要原预选V加密核查；不自动新增场景寻找优势。'))

page('八 可用于修订稿的英文结果段',
 *[p(x) for x in english.split('\n\n')[1:4]],
 p('上述文字与完整S4推导分开保存。英文回复草稿中的时态应在真正合入整稿后使用；目前交付的是可插入材料，未直接改写投稿原稿。'))

page('九 复现与验收记录',
 tbl(['产物','用途'],[
 ['inputs/protocol.json','结果前冻结的误差尺度、0.70°条件、半宽上限、三层检查和阴性结论规则'],
 ['inputs与sources','7条原始记录路径/哈希、315项冻结源码、完整配置与mpmath依赖身份'],
 ['results/*_export.mat','原MATLAB函数重建的数据、36维状态、参考/反馈/执行器裕度'],
 ['results/*_transfer.csv 与 *_bounds.csv','每个保存点的半宽、认证标记、覆盖及外向区间界'],
 ['summary_table 与 margin_summary','案例表、各有单位裕度的最小值和负值计数'],
 ['joint_failure_attribution 与 segments','未覆盖原因、保存点区段与最小半宽所在点'],
 ['S4_quantitative_transfer_draft.tex','完整有序推导与范围限定，供补充材料整合'],
 ['manuscript_insert_EN.txt','英文正文、结果与审稿回复草稿'],
 ['checks/independent_verification.json','70位独立核查、输入和源码完整性结果'],
 ],[2.75,3.95]),
 p('独立核查按每条轨迹12个等间隔索引、每层8个角点加中心，共2268次高精度状态评估；另检验6804个攻角导数值。它用于核对实现，全部误差盒的充分性来自区间包围，不来自这些离散角点。'),
 p('全样本MATLAB与独立Python相应量的重建数值差均小于8×10⁻¹³；所有逐点正半宽的区间通过标志已复核。新增轨迹积分为0，原7条轨迹和315项来源/冻结源码哈希全部保持一致。'),
 p('首次导出遇到旧Transparent缺少可选Y字段，修正后对可用字段核对并明确记录缺失标记；未生成或修改任何轨迹。材料包包含复现脚本、参数、重建状态与所有结果；原始大型MAT以路径和哈希引用，原论文工作副本不在本次重新发布。'))

# Minimal structured native OMML authoring; verified in the rendered Word export.
def run(s):
 z=OxmlElement('m:r');t=OxmlElement('m:t');t.text=s;z.append(t);return z
def many(node,terms):
 for x in terms:node.append(run(x) if isinstance(x,str) else deepcopy(x))
 return node
def subscript(a,b):
 z=OxmlElement('m:sSub')
 for tag,v in [('m:e',a),('m:sub',b)]:e=OxmlElement(tag);many(e,[v] if not isinstance(v,list) else v);z.append(e)
 return z
def sup(a,b):
 z=OxmlElement('m:sSup')
 for tag,v in [('m:e',a),('m:sup',b)]:e=OxmlElement(tag);many(e,[v] if not isinstance(v,list) else v);z.append(e)
 return z
def frac(a,b):
 z=OxmlElement('m:f')
 for tag,v in [('m:num',a),('m:den',b)]:e=OxmlElement(tag);many(e,[v] if not isinstance(v,list) else v);z.append(e)
 return z
def bar(a,pos='top'):
 z=OxmlElement('m:bar');pr=OxmlElement('m:barPr');v=OxmlElement('m:pos');v.set(qn('m:val'),pos);pr.append(v);z.append(pr);e=OxmlElement('m:e');many(e,[a]);z.append(e);return z
def hat(a):
 z=OxmlElement('m:acc');pr=OxmlElement('m:accPr');v=OxmlElement('m:chr');v.set(qn('m:val'),'̂');pr.append(v);z.append(pr);e=OxmlElement('m:e');many(e,[a]);z.append(e);return z
def formulas(name):
 s=subscript;ev=s('e','V');eg=s('e','γ');ec=s('e','χ');vf=s('V','f');vstar=s('V','*');af=s(hat('a'),'f');ac=s('a','c');mu=bar('μ','bot');lo=bar('a','bot');hi=bar('a');th=bar('Θ')
 if name=='box':return [[s('ℰ','ε'),' = {e : |',ev,'| ≤ ε,   |',eg,'|, |',ec,'| ≤ ',frac('ε',vstar),'}'],[vstar,' = 38.1 m/s']]
 if name=='attitude':return [['Θ = ',s('d','SO(3)'),'(',s('C','c'),',',s('C','r'),') + ',s('d','SO(3)'),'(',s('C','r'),',C) ≤ 0.70°']]
 if name=='margins':return [['ΔM(a) = K(2',vf,ev,' + ',sup(ev,'2'),')[',s('C','L'),'(a) + ',s('C','D'),'(a) tan a] + ΔX tan a'],
 [s('b','−,c'),' = ',s('b','−,f'),' + ΔN − ΔM(',s('a','−'),')'],[s('b','+,c'),' = ',s('b','+,f'),' + ΔM(',s('a','+'),') − ΔN']]
 if name=='root':return [[mu,' = inf ',s('M','a'),'(a; V, ',s('X','c'),') > 0'],['|',ac,' − ',af,'| ≤ ',frac(['sup |M(',af,'; V, ',s('X','c'),') − ',s('N','c'),'|'],mu),',     ',ac,' ∈ [',lo,',',hi,']']]
 if name=='thrust':return [['cos(',hi,' + ',th,') ≤ κ ≤ cos(max{0, ',lo,' − ',th,'})'],
 [s('T','raw'),' = ',frac([s('D','0'),'(V,',ac,') + ',s('X','c')],'κ'),',     ',bar(s('T','raw'),'bot'),' > 0,   ',bar(s('T','raw')),' < 60 N']]
 raise ValueError(name)

doc=Document();sec=doc.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11)
sec.top_margin=sec.bottom_margin=Inches(.65);sec.left_margin=sec.right_margin=Inches(.72);sec.footer_distance=Inches(.25)
for st,sz,bold in [('Normal',11,False),('Title',20,True),('Heading 1',15,True),('Heading 2',11.5,True)]:
 a=doc.styles[st];a.font.name='Microsoft YaHei';a.font.size=Pt(sz);a.font.bold=bold;a.font.color.rgb=RGBColor(0,0,0)
 a.paragraph_format.line_spacing=1.16;a.paragraph_format.space_after=Pt(7)
 rf=a.element.get_or_add_rPr().find(qn('w:rFonts'))
 if rf is None:rf=OxmlElement('w:rFonts');a.element.get_or_add_rPr().insert(0,rf)
 for x in ['ascii','hAnsi','eastAsia']:rf.set(qn('w:'+x),'Microsoft YaHei')
 if st!='Normal':a.paragraph_format.keep_with_next=True;a.paragraph_format.space_before=Pt(7)
for b in list(doc.styles.element.xpath('.//w:pBdr')):b.getparent().remove(b)
foot=sec.footer.paragraphs[0];foot.alignment=WD_ALIGN_PARAGRAPH.CENTER;fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');foot._p.append(fld)
plain=[]
for i,(title,blocks) in enumerate(P):
 if i:doc.add_page_break()
 doc.add_paragraph(title,'Title' if i==0 else 'Heading 1');plain.extend([title,''])
 for b in blocks:
  if b[0] in ['p','sub']:
   q=doc.add_paragraph(b[1],'Normal' if b[0]=='p' else 'Heading 2');plain.extend([b[1],''])
   if i==8:
    for r in q.runs:r.font.name='Calibri';r.font.size=Pt(11.5)
   continue
  if b[0]=='eq':
   for terms in formulas(b[1]):
    q=doc.add_paragraph();q.alignment=WD_ALIGN_PARAGRAPH.CENTER;q.paragraph_format.space_after=Pt(5)
    om=OxmlElement('m:oMath');many(om,terms);q._p.append(om)
   plain.append('[原生Word公式 '+b[1]+'；完整LaTeX推导见S4_quantitative_transfer_draft.tex]');continue
  if b[0]=='figure':
   q=doc.add_paragraph();q.alignment=WD_ALIGN_PARAGRAPH.CENTER;r=q.add_run();r.add_picture(str(H/'figures'/b[1]),width=Inches(b[2]));
   dp=q._p.xpath('.//wp:docPr')[0];dp.set('descr','原时变风Scalar逐保存点充分半宽、误差与姿态条件图');plain.append('[图 '+b[1]+']');continue
  _,heads,rows,widths=b;t=doc.add_table(rows=1,cols=len(heads));t.autofit=False;t.alignment=WD_TABLE_ALIGNMENT.CENTER
  for j,w in enumerate(widths):t.columns[j].width=Inches(w)
  for j,h in enumerate(heads):t.rows[0].cells[j].text=h
  for row in rows:
   cells=t.add_row().cells
   for j,v in enumerate(row):cells[j].text=v
  plain.append(' | '.join(heads));plain.extend(' | '.join(row) for row in rows);plain.append('')
  for ri,row in enumerate(t.rows):
   pr=row._tr.get_or_add_trPr();pr.append(OxmlElement('w:cantSplit'))
   if ri==0:pr.append(OxmlElement('w:tblHeader'))
   for ci,c in enumerate(row.cells):
    c.width=Inches(widths[ci]);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=c._tc.get_or_add_tcPr();mar=OxmlElement('w:tcMar')
    for side,val in [('top',80),('bottom',80),('left',95),('right',95)]:
     e=OxmlElement('w:'+side);e.set(qn('w:w'),str(val));e.set(qn('w:type'),'dxa');mar.append(e)
    pr.append(mar)
    if ri==0:e=OxmlElement('w:shd');e.set(qn('w:fill'),'E7EDF3');pr.append(e)
    for q in c.paragraphs:
     q.paragraph_format.line_spacing=1.1;q.paragraph_format.space_after=Pt(1)
     if i==4 and ci>1:q.alignment=WD_ALIGN_PARAGRAPH.CENTER
     for r in q.runs:r.font.size=Pt(9.5);r.font.bold=ri==0
  bord=OxmlElement('w:tblBorders')
  for edge in ['top','bottom','left','right','insideH','insideV']:
   e=OxmlElement('w:'+edge);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'D9D9D9');bord.append(e)
  t._tbl.tblPr.append(bord);doc.add_paragraph().paragraph_format.space_after=Pt(0)
doc.core_properties.title='Aerospace 原局部案例反馈裕度传递核查';doc.core_properties.author=''
doc.save(H/'Aerospace_WP3反馈裕度传递核查_v1_0.docx')
(H/'Aerospace_WP3反馈裕度传递核查_v1_0.txt').write_text('\n'.join(plain),encoding='utf-8')
print(json.dumps(dict(planned_pages=len(P),native_math_equations=len(doc.element.xpath('.//m:oMath')),tables=len(doc.tables)),ensure_ascii=False))
