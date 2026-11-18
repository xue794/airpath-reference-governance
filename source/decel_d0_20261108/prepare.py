"""Build D0 from the frozen causal six-DOF runner, locking all experiment inputs."""
from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent; R=H.parent
src=R/'contribution_probe_20261106/full6dof'; core=R/'revision_probe_20261106/source_snapshot/core'
reg=R/'decel_region_20261108'
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rep(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)
def put(n,s):(H/n).write_text(s,encoding='utf8')
def main():
    for k,v in json.loads((reg/'release_manifest.json').read_bytes()).items():assert h(reg/k)==v
    assert json.loads((reg/'results/decision.json').read_bytes())['passed']
    (H/'runs').mkdir(exist_ok=True);(H/'logs').mkdir(exist_ok=True)
    runner=(src/'run_cp6_warm.m').read_text()
    runner=runner[runner.index('clear cp6_governor_warm;'):]
    runner=rep(runner,"workspace=fileparts(fileparts(here));","workspace=fileparts(here);")
    runner=rep(runner,"addpath(fullfile(fileparts(here),'aoa'));","addpath(fullfile(workspace,'contribution_probe_20261106','aoa'));addpath(fullfile(workspace,'contribution_probe_20261106','full6dof'));addpath(fullfile(workspace,'decel_region_20261108'));")
    runner=rep(runner,'outdir=fullfile(here,tag);',"outdir=fullfile(here,'runs',tag);")
    runner=rep(runner,"'SpeedRise',8","'WindFun',@(tt)zeros(6,1),'SpeedRise',8")
    runner=runner.replace('cp6_snapshot_warm','d0_snapshot').replace('cp6_raw_margins','d0_raw_margins').replace('cp6_governor_warm','d0_governor')
    runner=runner.replace('d0_raw_margins(0,x0,C)','d0_raw_margins(0,x0,C,O)').replace('d0_raw_margins(tt,xx,C)','d0_raw_margins(tt,xx,C,O)')
    runner=rep(runner,"U=cl6_plant_input(xx(33:35),xx(36),'hold',C);","U=struct('delta',xx(33:35),'T',xx(36),'wind_fun',O.WindFun);")
    runner=rep(runner,'Results=cell(numel(methods),1);','d0_preflight(C,O,x0,method);\nResults=cell(numel(methods),1);')
    put('run_d0.m',"""function run_d0(method,maxStep)
assert(ismember(method,{'O','F','D','A','T','N10','Am'}) && ismember(maxStep,[.005,.01]));
Vlow=31;courseDegrees=0;turnRise=3;methods={method};
tag=sprintf('%s_%03d',method,round(1000*maxStep));
"""+runner)
    snap=(src/'cp6_snapshot_warm.m').read_text().replace('cp6_snapshot_warm','d0_snapshot').replace("cl6_airpath_state(t,xp,'hold',C)",'d0_air(t,xp,O.WindFun)').replace('cp6_governor_warm(xf,u0,C.R,method,O)','d0_governor(t,xf,u0,C.R,method,O)')
    put('d0_snapshot.m',snap)
    air=(core/'cl6_airpath_state.m').read_text()
    air=air[air.index('wI=z(1:3);'):]
    put('d0_air.m','function A=d0_air(t,xp,windFun)\nz=windFun(t);\n'+air)
    margins=(src/'cp6_raw_margins.m').read_text().replace('cp6_raw_margins(t,x,C)','d0_raw_margins(t,x,C,O)').replace("cl6_airpath_state(t,x(1:13),'hold',C)",'d0_air(t,x(1:13),O.WindFun)')
    put('d0_raw_margins.m',margins)
    gov=(src/'cp6_governor_warm.m').read_text()
    gov=gov[gov.index('% First try the updated'):];gov=rep(gov,"'rows',ix","'rows',ix,'cc',cc,'GG',GG,'backup',ub")
    put('d0_governor.m',"""function G=d0_governor(t,x,u0,R,method,O)
persistent previousSupport previousMethod
if any(strcmp(method,{'F','D','A'}))
    [C,ub]=moving_backup_affine(t,x,R);
else
    C=cp6_affine(x,R);ub=rg_backup_control(x,R);
end
switch method
 case {'O','F','A'}, ix=1:18;
 case 'D', ix=[1:13,16:18];
 case {'N10','Am'}, ix=1:17;
 case 'T', ix=[1:13,16:17];
 otherwise, error('Unexpected D0 arm');
end
cc=C.c(ix); GG=C.G(ix,:);
if any(strcmp(method,{'A','Am'}))
    Ac=cp_aoa_rows(x,R,'Am',C.alpha);
    cc(14:15)=Ac.c;GG(14:15,:)=Ac.G;
end
A=-GG;b=cc;Q=O.QP;Q.FeasibleStart=ub;
"""+gov)
    paths=[*H.glob('*.m'),*src.glob('*.m'),*core.glob('*.m'),*(R/'revision_probe_20261106/source_snapshot/W2D').glob('*.m'),reg/'release_manifest.json',reg/'moving_backup_affine.m',R/'contribution_probe_20261106/aoa/cp_aoa_rows.m',next(R.glob('Aerospace_有限时长*方案*docx'))]
    for mode in ('T','N10','Am'):
        for fine in (False,True):
            folder=(src/('V31_C60_fine' if fine else 'V31_C60')) if mode=='Am' else R/'endpoint_probe_20261106/full6dof'/('V31_C60_warm_fine' if fine else 'V31_C60_warm')
            paths.append(folder/(mode+'.mat'))
    config=dict(scope='D0 only',speed_start=3,ramp=8,target=31,turn=0,horizon=35,wind=0,noise=0,steps=[.01,.005],sample=.05,arms=['O','F','D','A','T','N10','Am'],windows=[[3,16],[14,20],[25,35]],tracking_min=.606,force_min_N=1,drift_multiplier=5,speed_floor=.0001,force_floor=.0000001,same_state_scaled_tol=1e-7,intervention_scaled_tol=1e-7,prefix_tol_state_abs=1e-5,inputs={p.relative_to(R).as_posix():h(p) for p in paths})
    lock=H/'protocol.json'
    if lock.exists():assert json.loads(lock.read_bytes())==config
    else:lock.write_text(json.dumps(config,indent=2),encoding='utf8')
if __name__=='__main__':main()
