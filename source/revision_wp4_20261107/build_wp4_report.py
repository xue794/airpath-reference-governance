from pathlib import Path
import json,math
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
H=Path(__file__).resolve().parent
S=json.loads((H/'results/summary_table.json').read_text(encoding='utf-8'))
by={x['method']+'_'+x['case']:x for x in S}
case={'hold':'H 静风','varying':'V 时变风','noisy':'N 测量误差'}
name={'Transparent':'Transparent','Proposed':'Scalar','MatchedQP':'Matched-QP'}
P=[]
def p(x):return ('p',x)
def sub(x):return ('sub',x)
def tab(h,r,w):return ('table',h,r,w)
def fig(x,w=6.65):return ('figure',x,w)
def page(t,*blocks):P.append((t,blocks))
def down(x):return f'{math.floor(x*1e6)/1e6:.6f}' if x>0 else '0 未认证'

english='''Results paragraph for integration into the revised manuscript

The original local comparison now includes all nine case-method cells: hold, varying wind, and varying wind with the original deterministic measurement-error signals, each evaluated with Transparent, Scalar, and Matched-QP. Seven trajectories were reused and only the two missing Transparent trajectories were integrated. The added runs retain the frozen model, initialization, gains, wind and reference profiles, inversion interval, and numerical settings. Transparent applies the original nominal reference input and retains the same downstream controller and actuator model. All nine trajectories complete their prescribed horizons, with 18,009 retained samples in total.

In hold, the three methods agree, with a position RMS of 0.027789 m. In the noisy case, the position RMS values are 2.416064, 7.750469, and 8.409029 m for Transparent, Scalar, and Matched-QP, respectively. Transparent violates the 3 N reference-thrust reserve at 157 of 2,401 saved points and leaves the selected backup set at 2,339 points. Its reference thrust is negative at 17 saved points from 26.60 to 27.40 s, with a minimum of -0.277596 N. These are sampled reference-layer violations: the raw feedback thrust remains between 10.681875 and 27.817714 N, and the reconstructed feedback inversion margins, modeled actuator limits, and checked study-domain limits are satisfied at every saved point. No saved-point actuator limiting or active allocation bound is detected. Scalar and Matched-QP preserve the audited reference constraints on these trajectories. The comparison therefore documents the cost of reference admissibility and selected-set containment; it does not demonstrate that the observed tracking loss is necessary for physical protection or that endpoint constraints provide a benefit beyond an existing thrust constraint.

Applying the unchanged WP3 interval definitions gives joint-condition inclusion of 1,201/1,201 for Transparent hold and 1,567/2,401 for Transparent noisy. In the latter, the 17 negative-reference-thrust points have no positive joint certificate; another 817 points have a positive sufficient half-width that does not include the observed air-path error. Endpoint and geometry checks include the observed errors at all retained points, and the original 0.70-degree attitude-sum condition is observed throughout. These finite-record statements neither establish error-box invariance nor certify all between-sample times. The original seven-cell numerical results are unchanged.

Response to Reviewer 1 after the additions are merged

The revised comparison separates reference design reserves, raw inversion feasibility, backup-set membership, feedback commands, and modeled actuator states. In the original varying-wind Transparent record, 155/2,401 points violate the reference-thrust reserve and 2,339/2,401 lie outside the selected backup set; no retained-point feedback or modeled physical-limit violation is found. The added noisy record also contains 17 negative-reference-thrust samples, while its raw feedback thrust remains feasible. We therefore describe the tracking tradeoff at the reference layer and have removed any inference that these examples alone demonstrate necessary physical protection.

Response to Reviewer 3 after the additions are merged

This audit supports a reproducible, architecture-specific reference admissibility comparison and quantifies its local feedback-transfer limitations. It does not isolate an endpoint-specific advantage over a matched angle-of-attack or thrust-limited construction. Accordingly, closed-form selection, absence of online optimization, and recursive feasibility are not presented as standalone novelty claims, and the present local examples are not used as proof of such an independent engineering advantage.
'''
(H/'manuscript_insert_EN.txt').write_text(english,encoding='utf-8')

page('Transparent 记录补齐与分层约束核查',
 p('Aerospace 修订 WP4 完成报告  2026年11月7日'),
 p('原方案中的两条记录缺口已经补齐。Transparent hold 完成60 s，Transparent noisy完成120 s；两者均无异常退出。与七条既有轨迹合并后，原九单元矩阵共有18,009个保存点。'),
 p('新增证据确认了参考层的约束代价。Transparent noisy跟踪更准，但有157个点失去3 N推力预留，其中17个点的参考推力为负；2339个点在所选备份集外。实际反馈推力和反演根仍可行，因此这些结果不能直接证明牺牲跟踪精度换来了必要的物理保护。'),
 tab(['原案例','Transparent\n位置RMS m','Scalar\n位置RMS m','Matched-QP\n位置RMS m'],[
 [case[c]]+[f"{by[m+'_'+c]['RMS_grid_m']:.6f}" for m in ['Transparent','Proposed','MatchedQP']] for c in ['hold','varying','noisy']
 ],[1.6,1.65,1.65,1.8]),
 p('RMS采用原0.05 s保存网格的全程位置误差定义。H为60 s，V/N为120 s；N使用原确定性测量误差波形，不是随机重复试验。数据另附时间积分RMS、后段RMS和峰值。暂无独立任务误差指标，不给出任务合格或不合格结论。'),
 sub('本次完成范围'),
 tab(['项目','执行结果'],[
 ['补跑与复用','仅新增Transparent H/N两次积分；复用其余七条，不重跑Monte Carlo。'],
 ['保持条件','原十八行、备份集合参数0.10、原矩阵与增益、原0–9°反演区间、10 N与3 N预留。'],
 ['公平性检查','相同初态、模型与下游控制；仅Transparent参考输入使用原名义值。积分前代数核对通过。'],
 ['研究边界','未启用新的恢复备份、认证集合扩展、可选L场景或V加密积分。'],
 ],[1.3,5.4]),
 p('对论文贡献的含义：现有专用构造可以作为可计算的参考可行性接口描述；这批结果仍未分离出端点约束相对已有推力限制或强匹配攻角约束的独立收益。'))

page('一 九单元的参考约束记录',
 p('下表逐层保留参考状态的越界情况。3 N预留不满足、参考推力低于0 N、离开备份集分别计数。负参考推力是所选虚拟参考所要求的推力，不是飞机实际施加了负推力。'),
 tab(['案例','方法','参考推力\n最小 N','低于3 N\n点数','低于0 N\n点数','备份集外\n点数'],[
 [case[x['case']],name[x['method']],f"{x['reference_T_min_N']:.6f}",str(x['reference_thrust_reserve_violations']),str(x['reference_other_state_violations']),str(x['backup_outside_samples'])] for x in S
 ],[1.05,1.35,1.25,1,1,1.05]),
 p('本批记录中，参考层除推力预留、N案例的0 N下界和备份集外，没有其他状态约束越界。上表“低于0 N”计数已与逐行裕度表交叉确认；一般情况下不能把其他参考状态违约都归入该列。'),
 sub('导数条件与状态条件分开'),
 p('Transparent noisy的非备份导数条件有163个点不满足，均对应推力下界行；备份导数条件有1810个点不满足。其首次保存时刻分别为23.75 s与2.05 s。当前状态实际首次失去3 N预留为24.80 s，首次在备份集外为3.10 s。导数条件为负不等于同一时刻状态条件已经为负。'),
 sub('反演与模型物理量'),
 p('九条记录均未发现保存点上的参考或反馈反演失败。Transparent noisy的参考原始端点最小裕度为208.60 N和96.26 N，反馈端点为205.38 N和96.11 N；即使扣除10 N也保持正值。其参考根约4.242–5.870°，反馈命令根约4.166–5.824°，真实攻角约4.167–5.829°，三者分开记录。'))

page('二 新发现的参考负推力',
 p('Transparent noisy在26.60–27.40 s的17个保存点上要求负参考推力，最小值为−0.277596 N，出现在27.00 s。这些首末保存时刻不是经过事件定位的连续时间边界。'),
 fig('noisy_negative_reference.png'),
 p('图1  左图将原记录的参考推力下界局部放大，红点为17个负值保存点。右图使用相同时间范围显示限幅前原始反馈推力和实际执行器状态。连线只帮助读取采样记录。'),
 tab(['量','Transparent noisy全程结果'],[
 ['参考推力','−0.277596至38.776130 N；157点低于3 N，其中17点低于0 N。'],
 ['原始反馈推力','10.681875至27.817714 N；无0–60 N或3–57 N保存点越界。'],
 ['反演适用域','所有保存点的原0–9°参考和反馈根均存在，端点保持正裕度。'],
 ['备份集合','2339/2401点在集合外；备份函数值与集合参数之比最大约20.316。'],
 ],[1.35,5.35]),
 p('70位独立计算重新求解全部17个点及前后相邻两点的参考根和推力，最大推力差约1.03×10⁻¹⁴ N，负值分类一致。没有扩大求根区间、补偿异常或为获得有利结果重调波形。'),
 p('九条记录的原始反馈推力、已检查反馈几何裕度、模型执行器幅值与速率，以及研究采用的实际空速、航迹角与攻角界限均无保存点越界。执行器限幅、速率限制和分配器边界激活计数均为0。研究界限不是经过外部认证的飞机硬件包线。'),
 p('反馈控制使用实际空气状态及跟踪误差修正，故反馈命令不必等于虚拟参考所需推力。本例中，参考不可行性并未演变为记录内的反馈或执行器违约。'))

page('三 沿用 WP3 的局部传递核查',
 p('数值定义保持不变：角度误差乘以38.1 m/s，与空速误差组成逐分量误差盒；姿态条件为原S4两段姿态误差之和不超过0.70°；搜索上限1 m/s，二分20次，40位区间运算。端点及几何条件先成立，随后包围隐式根和原始反馈推力。'),
 tab(['案例','方法','J最小半宽\nm/s','正半宽\n点数','联合包含\n点数','记录比例'],[
 [case[x['case']],name[x['method']],down(x['J_min_mps']),str(x['J_positive']),f"{x['J_inclusion']}/{x['samples']}",f"{x['J_percent']:.2f}%"] for x in S
 ],[1.05,1.35,1.2,.85,1.25,1.0]),
 p('正最小半宽向下保留6位小数；精确通过值保存在逐点CSV。0表示未获得正证书，不是“认证半径为零”。联合包含要求实际气路误差进入该保存参考点的充分盒，同时满足姿态条件。比例是记录描述，不是安全概率或连续时间覆盖率。'),
 p('九条记录的端点E层和几何G层均达到1 m/s搜索上限，并包含全部保存点的实际误差；所有点都满足记录中的0.70°姿态和条件。这个搜索上限不是最大允许半宽，也没有证明误差未来保持在盒中。'),
 p('Transparent noisy的17个J层零返回值与负参考推力保存点完全一致，原点即未通过联合判据。另有817点虽有正半宽，实际气路误差却超出所算盒；合计834个未包含点。对实际误差盒的推力充分包围检查中，803点未通过下界，31点未通过上界。不能由包围失败推出实际推力越界。'),
 p('WP3原七条记录的数字完全保持不变。新增hold为1201/1201，新增noisy为1567/2401。以事后观测的直接姿态角替换固定条件时，noisy会得到1568点，但这只作诊断，主结果保持1567；不能用事后选择替代先验保证。'))

page('四 这些结果应怎样进入大修稿',
 tab(['审稿问题','可写入的证据','解释边界'],[
 ['Transparent违反哪些约束','V/N的参考推力预留和备份集计数；N另有17点负参考推力','不写成实际执行器越界、飞机失稳或反馈求根失败。'],
 ['跟踪损失换来了什么','Scalar与Matched-QP保持已审计参考约束及所选集合；Transparent的V/N位置RMS更小','现有例子未证明这些损失为物理保护所必需。'],
 ['局部理论范围是否可核查','原九单元的E/G/J层数值与包含计数，公开17个无正J证书点','不宣称整个备份集统一半径、误差管不变或采样间保证。'],
 ['专用端点是否有独立收益','原局部轨迹的端点裕度较充分；主要缺口在推力与集合条件','新增负参考推力由推力约束即可直接识别，不能归功于端点表示独有优势。'],
 ],[1.15,2.8,2.75]),
 sub('正文与补充材料的有限修改'),
 p('结果节用九单元矩阵替换Transparent缺项，并增加分层违约表及图1。S9删除复用旧配置全约束违约比例的入口，改为本次统一配置的逐行统计。S4在WP3定义不变的基础上补两条记录，明确Transparent noisy的17个未认证点。'),
 p('讨论段保留计算与匹配QP对照的有限结论，将贡献表述集中于耦合接口的可计算约束和条件性传递核查。“闭式更新”“无在线优化”“递归可行”继续作为实现与理论特征，不能单独承担新颖性。'),
 p('本次补缺完成后，可以进入WP5整合稿件。仅凭本次数据不能保证审稿人3接受核心贡献，但可以准确回答Transparent究竟违反了什么，并去掉现有证据不支持的工程收益表述。研究主线不转向新的恢复备份设计。'),
 sub('原方案中的可选项'),
 p('本次冻结清单未启用L场景或V加密积分。现有结论严格限定为保存点记录，负参考推力已经独立高精度核查，不将事件精确定位或连续时间保证纳入本轮。原九单元缺口关闭，可选分支不自动变成新增研究任务。'))

page('五 可插入修订稿的英文结果',
 *[p(x) for x in english.split('\n\n')[1:4]],
 p('英文回复审稿人1与3的建议段另存于manuscript_insert_EN.txt。回复中“已加入”“已删除”等表述应在真正合并正文与补充材料后采用；本阶段保留独立插入材料。'))

page('六 复现记录与阶段验收',
 tab(['核查项','结果与定位'],[
 ['归档身份','检查65个ZIP目录及相关文件。近似命名BASELINE_020 H/N属于W1C旧Scalar，不能用作Transparent。详见archive_inventory.json。'],
 ['运行前冻结','protocol_manifest.json和run_fingerprint.json绑定参数、源文件、驱动和两个冻结输入。318个源文件快照。'],
 ['积分设置','MATLAB R2018a；ode45；RelTol 1e−8，AbsTol 1e−10，MaxStep 0.01 s；保存间隔0.05 s。'],
 ['新增记录','runs目录保存Transparent hold/noisy原始36状态、输出、状态清单及过程日志。两次均完整结束，无异常。'],
 ['hold一致性','Transparent与原Scalar完整36状态记录最大绝对差1.53×10⁻¹³；满足原1e−7比较容差。'],
 ['统一分层重建','原七条RMS、推力极值和联合计数不变。旧Transparent V存档的集合参数0.20仅在审计时统一为0.10；其输入不使用该参数。新H/N固定0.10。'],
 ['区间独立核对','新增两条的648个误差状态、1944个导数值及216个根通过70位核查；并检查全部3602点保存证书。'],
 ['完整重建差','合并九条记录的最大独立重建差小于8×10⁻¹³；负参考推力另作19点针对性高精度核对。'],
 ],[1.25,5.45]),
 p('CSV中的负裕度判定容差为1e−8；导数行保留原始值及归一化诊断。数学充分证书仍使用严格正条件，未以这一记录容差放宽。参考、反馈、执行器及研究界限分开输出，单位在margin_summary中注明。'),
 p('技术修正：九条目录的JSON条目字段不完全一致，MATLAB将其解码为cell。首轮统一审计在读取第一条前停止；修正只涉及目录解析，随后全量重建成功。轨迹、冻结输入和源控制律没有改变，也没有重做积分。失败日志保留。'),
 p('阶段状态：WP4原矩阵补缺已完成；WP5整合正文和补充材料待执行。报告、可插入文本、逐点记录、核查脚本与复现说明一并归档。既有投稿文件、第一阶段冻结基线及WP3成果均保持原样。'))

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
    heading=doc.add_paragraph(title,'Title' if i==0 else 'Heading 1')
    if i:heading.paragraph_format.page_break_before=True
    plain.extend([title,''])
    for b in blocks:
        if b[0] in ['p','sub']:
            q=doc.add_paragraph(b[1],'Normal' if b[0]=='p' else 'Heading 2');plain.extend([b[1],''])
            if i==5 and not b[1].startswith('英文回复'):
                for r in q.runs:r.font.name='Calibri';r.font.size=Pt(11.5)
            continue
        if b[0]=='figure':
            q=doc.add_paragraph();q.alignment=WD_ALIGN_PARAGRAPH.CENTER;q.add_run().add_picture(str(H/'figures'/b[1]),width=Inches(b[2]));q.paragraph_format.keep_with_next=True
            q._p.xpath('.//wp:docPr')[0].set('descr','Transparent noisy参考负推力保存点及对应的反馈与执行器推力');plain.append('[图 '+b[1]+']');continue
        _,heads,rows,widths=b;t=doc.add_table(rows=1,cols=len(heads));t.autofit=False;t.alignment=WD_TABLE_ALIGNMENT.CENTER
        for j,w in enumerate(widths):t.columns[j].width=Inches(w)
        for j,h in enumerate(heads):t.rows[0].cells[j].text=h
        for row in rows:
            cells=t.add_row().cells
            for j,v in enumerate(row):cells[j].text=str(v)
        plain.append(' | '.join(heads));plain.extend(' | '.join(map(str,row)) for row in rows);plain.append('')
        for ri,row in enumerate(t.rows):
            pr=row._tr.get_or_add_trPr();pr.append(OxmlElement('w:cantSplit'))
            if ri==0:pr.append(OxmlElement('w:tblHeader'))
            for ci,c in enumerate(row.cells):
                c.width=Inches(widths[ci]);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=c._tc.get_or_add_tcPr();mar=OxmlElement('w:tcMar')
                for side,val in [('top',75),('bottom',75),('left',90),('right',90)]:
                    e=OxmlElement('w:'+side);e.set(qn('w:w'),str(val));e.set(qn('w:type'),'dxa');mar.append(e)
                pr.append(mar)
                if ri==0:e=OxmlElement('w:shd');e.set(qn('w:fill'),'E7EDF3');pr.append(e)
                for q in c.paragraphs:
                    q.paragraph_format.line_spacing=1.1;q.paragraph_format.space_after=Pt(1)
                    for r in q.runs:r.font.size=Pt(9.5);r.font.bold=ri==0
        border=OxmlElement('w:tblBorders')
        for edge in ['top','bottom','left','right','insideH','insideV']:
            e=OxmlElement('w:'+edge);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'D9D9D9');border.append(e)
        t._tbl.tblPr.append(border);doc.add_paragraph().paragraph_format.space_after=Pt(0)
doc.core_properties.title='Transparent 记录补齐与分层约束核查';doc.core_properties.author=''
doc.save(H/'Aerospace_WP4记录补齐与分层核查_v1_0.docx')
(H/'Aerospace_WP4记录补齐与分层核查_v1_0.txt').write_text('\n'.join(plain),encoding='utf-8')
print(json.dumps(dict(planned_pages=len(P),tables=len(doc.tables)),ensure_ascii=False))
