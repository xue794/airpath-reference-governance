from pathlib import Path
import json, hashlib
h=Path(__file__).resolve().parent
root=h.parent
old=root/'decel_d0_20261108'
def adapt(source, dest, changes):
    text=(old/source).read_text(encoding='utf8')
    for a,b in changes:
        assert text.count(a)==1,a
        text=text.replace(a,b)
    (h/dest).write_text(text,encoding='utf8')
adapt('run_d0.m','run_confirm.m',[
('function run_d0(method,maxStep)','function run_confirm(scenario,method,maxStep)'),
('Vlow=31;courseDegrees',"""assert(ismember(scenario,{'C1','C2','C3'}));
assert(strcmp(scenario,'C1') || ismember(method,{'O','F','D','A'}));
if strcmp(scenario,'C1'),Vlow=31;else,Vlow=31.5;end
if strcmp(scenario,'C3'),windFun=@confirm_wind;else,windFun=@(tt)zeros(6,1);end
courseDegrees"""),
("tag=sprintf('%s_%03d',method","tag=sprintf('%s_%s_%03d',scenario,method"),
('source=fullfile','addpath(fullfile(workspace,\'decel_d0_20261108\'));\nsource=fullfile'),
("'WindFun',@(tt)zeros(6,1),'SpeedRise',8","'WindFun',windFun,'SpeedRise',10"),
('d0_preflight(C,O,x0,method);','confirm_preflight(C,O,x0,method,scenario);')])
adapt('analysis/replay_d0.m','analysis/replay_confirm.m',[
('function replay_d0(method,step)','function replay_confirm(scenario,method,step)'),
('addpath(here,fullfile','addpath(here,fullfile(root,\'decel_d0_20261108\'),fullfile'),
("sprintf('%s_%03d',method","sprintf('%s_%s_%03d',scenario,method"),
("'FDKKT','FAKKT'}","'FDKKT','FAKKT','airSpeed','refAirSpeed','windN','windE','windD','windDotN','windDotE','windDotD','whatN','whatE','whatD'}"),
("abs(S.Cmd.aChi),dr]","abs(S.Cmd.aChi),dr,S.Air.V,xf(1),S.Air.wI',S.Air.wDotI',xx(17:19)']"),
('D0_REPLAY_OK','CONFIRM_REPLAY_OK')])
adapt('analysis/chart_events.m','analysis/confirm_chart_events.m',[
('function chart_events()','function confirm_chart_events()'),('addpath(here,study,','addpath(here,study,fullfile(root,\'decel_d0_20261108\'),'),
("fullfile(study,'d0_snapshot.m')","fullfile(root,'decel_d0_20261108','d0_snapshot.m')")])
adapt('analysis/free_ablation.m','analysis/confirm_free.m',[
('function free_ablation()','function confirm_free()'),('addpath(h,fullfile','addpath(h,fullfile(r,\'decel_d0_20261108\'),fullfile(r,\'decel_region_20261108\'),fullfile'),
("['N10_' id{1}]","['C1_N10_' id{1}]")])
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
cfg=json.loads((old/'protocol.json').read_bytes())
cfg.update(scope='C1 C2 C3 confirmation',cases={'C1':{'target':31,'wind_gain':0},'C2':{'target':31.5,'wind_gain':0},'C3':{'target':31.5,'wind_gain':0.1}},ramp=10,arms={'C1':['O','F','D','A','T','N10','Am'],'C2':['O','F','D','A'],'C3':['O','F','D','A']},wind='gain * smoothstep5(clip(t/3,0,1)) * wind_estimator_test_wind(t,varying), derivative by product rule',stage_gate='D0 released; optional Scalar only after strong endpoint and tracking outcome on all confirmation cases')
cfg.pop('target')
for p in [old/'protocol.json',old/'release_manifest.json',root/'tmp/decel_plan/content.py',next(root.glob('*有限时长减速机动研究方案*.docx')), *h.glob('*.m'),*[h/'analysis'/n for n in ['replay_confirm.m','confirm_chart_events.m','confirm_free.m']]]:
    cfg['inputs'][p.relative_to(root).as_posix()]=digest(p)
for path,sha in cfg['inputs'].items():assert digest(root/path)==sha,path
dest=h/'protocol.json'
value=json.dumps(cfg,ensure_ascii=False,indent=2)
if dest.exists():assert json.loads(dest.read_bytes())==cfg
else: dest.write_text(value,encoding='utf8')
