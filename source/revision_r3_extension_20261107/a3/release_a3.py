"""Render the A3 evidence and package the MATLAB/Python replay paths with verified inputs."""
import hashlib
import json
import csv
from html import escape
from zipfile import ZipFile, ZIP_DEFLATED

from prepare import HERE, ROOT, A1, A2, sha

def get_csv(name):
    with (HERE/name).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def table(head,data):
    def tags(kind,cells): return "".join(f"<{kind}>{escape(str(s))}</{kind}>" for s in cells)
    return '<div class="scroll"><table><thead><tr>'+tags("th",head)+'</tr></thead><tbody>'+''.join(
        '<tr>'+tags("td",r)+'</tr>' for r in data)+'</tbody></table></div>'

def format_num(row,key,places=4,sign=False):
    if row[key] in ("",None):return "—"
    return format(float(row[key]),f"{'+' if sign else ''}.{places}f")

def report():
    z=json.loads((HERE/"a3_verified.json").read_bytes())
    r=get_csv("run_summary.csv"); c=get_csv("contrasts.csv"); fine=[x for x in r if x["step"]=="fine"]
    f=format_num
    words={"higher":"N10 较高","lower":"N10 较低","unresolved":"未达分辨门槛","window_unavailable":"无共同窗口"}
    html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A3 闭环反馈端点裕度对照</title><style>
body{max-width:1160px;margin:0 auto;padding:30px min(4vw,52px) 64px;font:16px/1.7 "Microsoft YaHei",sans-serif;color:#172a33}
h1{font-size:30px;line-height:1.3}h2{font-size:21px;color:#195169;margin:30px 0 8px}p{margin:8px 0}
table{border-collapse:collapse;width:100%;font:14px/1.55 "Microsoft YaHei",sans-serif;font-variant-numeric:tabular-nums}
td,th{text-align:right;padding:6px 8px;border-bottom:1px solid #d5e1e2;white-space:nowrap}th{background:#e7f0f1}
td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}a{color:#13556c}@media print{.scroll{overflow:visible}tr{break-inside:avoid}}</style><main>
<h1>A3 闭环反馈上端点裕度对照</h1>'''
    html+=(HERE/"interpretation.html.inc").read_text(encoding="utf-8")
    html+='''<h2>同一初态的三种闭环</h2><p>按 A2 有确认点的样本编号升序取前三例。N10 使用原 18 行，D16 删除上下端点的两条导数行，Am 用等 10 N 预留的移动攻角导数行替换它们，保留边界导数。三组的 QP 目标、备份椭球与控制、3–57 N 参考推力限制、反馈控制器、36 维初态和各例的扰动一致。原始反馈推力受独立的 0–60 N 执行器指令限幅。</p>
<p>全部轨迹从 0 s 出发。N10 原步长取自 A1，逐点重建指令、推力及参考和反馈端点间隙。其余 15 条使用 MATLAB R2018a 的 ode45 重新积分。粗、细最大步长为 0.01、0.005 s，输出间隔均为 0.05 s，相对与绝对容差为 10⁻⁸、10⁻¹⁰。比较窗口含原 N10 首次 A2 确认点，止于同一编号三个配置两种步长共有的最后保存点。</p>'''
    html+='<h2>预定主指标</h2><p>每组先在窗口内求反馈上端点力间隙 g₊ = M(9°, V, X) − N 的最小保存值，再分别作差。正号表示 N10 的最小裕度较大。表中的门槛为 5 × max(两组各自的粗细漂移，10⁻⁷ N)。粗细差值须同号且两个绝对值均严格越过门槛。此门槛只判断计算分辨率。</p>'
    def sig(k,col,signed=False):
        return format(float(k[col]),("+" if signed else "")+".6g") if k[col] else "—"
    html+=table(["编号","对照","窗口 (s)","粗步差 (N)","细步差 (N)","N10 漂移 (N)","对照漂移 (N)","5 倍门槛 (N)","判定"],
        [[k["id"],k["contrast"],f(k,"start_s",2)+"–"+f(k,"common_end_s",2),sig(k,"coarse_delta_N",True),
          sig(k,"fine_delta_N",True),sig(k,"N10_drift_N"),sig(k,"other_drift_N"),sig(k,"five_x_gate_N"),words[k["classification"]]] for k in c])
    html+='<h2>逐例细步结果</h2><p>最小上端点和出现时刻按上述窗口计算。下表的位置 RMS、参考力间隙、推力范围及备份成员关系覆盖自 0 s 起的完整共同区间。RMS 取各保存点三维位置误差平方范数的均值再开方。</p>'
    html+=table(["编号","配置","原终点 (s)","最小上端点 (N)","发生时刻 (s)","位置 RMS (m)","参考下端最小 (N)","参考上端最小 (N)","反馈下端最小 (N)","最大 V_B/ρ_B"],
        [[k["id"],k["arm"],f(k,"original_end_s",2),f(k,"primary_min_upper_N",6),f(k,"primary_argmin_s",2),f(k,"full_common_RMS_m",3),
          f(k,"min_ref_raw_lower_N",3),f(k,"min_ref_raw_upper_N",3),f(k,"min_fb_raw_lower_N",3),f(k,"max_VB_over_rho",9)] for k in fine])
    def rng(k,col):return f(k,"min_"+col,3)+" 至 "+f(k,"max_"+col,3)
    html+=table(["编号","配置","参考推力 (N)","反馈原始请求 (N)","执行器实际推力 (N)","请求裁剪点","舵面位置最小余量 (rad)","舵速最小余量 (rad/s)","推力速率最小余量 (N/s)"],
        [[k["id"],k["arm"],rng(k,"ref_T_N"),rng(k,"T_raw_N"),rng(k,"T_act_N"),k["clip_samples_common"],f(k,"min_surface_margin_rad",5),
          f(k,"min_surface_rate_margin_radps",5),f(k,"min_T_rate_margin_Nps",4)] for k in fine])
    html+='<h2>移动攻角界的备份核查</h2><p>下表在 Am 实际经过的全部共同保存点，将同一个备份控制 u_B 代入下、上攻角导数行。角度状态余量以 rad 计，导数行余量以 rad/s 计。使用独立求根和沿备份向量场的中心差分核查移动边界导数。这项检查的范围为表列轨迹；原力端点全椭球证书的适用范围仍由原行定义。</p>'
    html+=table(["编号","最大步长 (s)","下界导数最小 (rad/s)","上界导数最小 (rad/s)","下界状态最小 (rad)","上界状态最小 (rad)","最大归一化备份违反量"],
        [[k["id"],"0.01" if k["step"]=="coarse" else "0.005",f(k,"min_Am_backup_lo_radps",6),f(k,"min_Am_backup_hi_radps",6),
          f(k,"min_Am_h_lo_rad",6),f(k,"min_Am_h_hi_rad",6),format(float(k["max_backup_violation"]),".3g")] for k in r if k["arm"]=="Am"])
    html+=f'''<p>完整原时域共 {z["points"]:,} 个配置保存点，其中 {z["completed"]}/{z["cases"]} 条配置记录到 120 s。表中所有极值均按 0.05 s 输出网格统计。完整数据按每条轨迹保存结束状态、全部 36 维状态、控制与力裕度序列及数值残差。</p>
<p><a href="contrasts.csv">粗细步长对照</a> · <a href="run_summary.csv">18 组指标与原终点</a> · <a href="model_limits.csv">模型界触及台账</a> · <a href="thrust_clipping.csv">裁剪区间</a> · <a href="a3_verified.json">独立核查误差</a> · <a href="README_复算说明.txt">复算说明</a></p></main></html>'''
    (HERE/"A3_闭环反馈裕度对照.html").write_text(html,encoding="utf-8")

def release():
    report()
    files={ROOT/p for p in json.loads((HERE/"input_hashes.json").read_bytes())}
    files.update([A1/"release_manifest.json",A2/"release_manifest.json",HERE/"input_hashes.json",HERE/"a3_verified.json",
        HERE/"batch.log",HERE/"interpretation.html.inc"])
    for pat in ["*.py","*.csv","*.html","*.txt","runs/*_states.mat","runs/*_samples.mat","runs/*_status.json","identity/*.csv"]:
        files.update(HERE.glob(pat))
    manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)}
    manifest_path=HERE/"release_manifest.json"
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    path=ROOT/"Aerospace_A3闭环反馈裕度对照材料包_v1_0.zip"
    temp=path.with_suffix(".partial.zip")
    with ZipFile(temp,"w",ZIP_DEFLATED,compresslevel=3) as archive:
        for p in sorted(files|{manifest_path}):
            archive.write(p,p.relative_to(ROOT).as_posix())
    with ZipFile(temp) as archive:
        assert len(archive.namelist())==len(files)+1
        for key,record in manifest.items():
            with archive.open(key) as f: assert hashlib.file_digest(f,"sha256").hexdigest()==record["sha256"]
    temp.replace(path)
    print(len(files)+1,"files verified",path.stat().st_size,"bytes")

if __name__=="__main__":
    release()
