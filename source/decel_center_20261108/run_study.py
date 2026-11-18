"""Lock the post-confirmation comparator, then run four independent MATLAB jobs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import hashlib, json, shutil, subprocess, sys
H=Path(__file__).resolve().parent
W=H.parent
C=W/'decel_confirm_20261108'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def freeze():
    from build_adapter import build
    build()
    old=json.loads((C/'protocol.json').read_bytes())
    for path,h in old['inputs'].items():assert sha(W/path)==h,path
    archive=json.loads((C/'release_manifest.json').read_bytes())
    files=[C/'protocol.json',C/'release_manifest.json',W/'decel_region_20261108/results/decision.json',*sorted(H.glob('*.m')),H/'build_adapter.py',H/'run_study.py']
    for p in [C/'results/run_metrics.csv',*(C/'runs').glob('C[12]_F_*/*'),*(C/'runs').glob('C[12]_O_*/replay.csv')]:
        assert archive[p.relative_to(C).as_posix()]==sha(p)
        files.append(p)
    spec=dict(status='Post hoc supplementary experiment after opening C1/C2; previously observed F and nominal center comparison.',cases={'C1':31,'C2':31.5},steps=[.01,.005],
        arm='Direct analytic center (38.1 to 31.25 m/s, starts 3 s, 22 s quintic), level and straight; actual six-DOF closed loop.',
        inputs_to_controller='Use exact center and analytic rates. Bypass governor and upstream command selection. Original guidance and observer evolve for diagnostics. Same saved F configuration, trim, physical model and downstream feedback.',
        environment='Both calm, no measurement noise. Ground requests start at 3 s with 10 s ramp, horizon 35 s; sample interval .05 s.',
        events='Downward first crossing of either original endpoint in reference or feedback terminates. Raw guidance speed and horizontal chart thresholds are passive logged events, as direct controller does not consume that command. For unequal lengths report common span only.',
        primary='Actual horizontal groundspeed time-averaged absolute error [3,16]; compare Direct minus F. Descriptive reference-only counterpart and signed two-difference decomposition using historical O over same complete window.',
        companion=['35 s sampled 3D position RMS and peak','feedback upper and lower minima at 14–20 s and until stop with times','mean actual and governed ground speed at 25–35 s, signed target bias and C2-C1 response','E/G/J sampled coverage and actuators; distinguish full-domain proof and observed samples'],
        interpretation='Same sign and difference >5*max(coarse/fine drift in either arm, floor), with floors 1e-4 m/s or 1e-7 N. Keep old research scales .606 m/s and 1 N. Position descriptive. No mission pass thresholds. Either direction and null outcomes retained.',
        hashes={p.relative_to(W).as_posix():sha(p) for p in files})
    p=H/'protocol.json'
    if p.exists():
        old=json.loads(p.read_bytes());old.pop('registered_local');assert old==spec
    else:
        spec['registered_local']=datetime.now().astimezone().isoformat();p.write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
def job(s):
    case,step=s;tag=f'{case}_{round(step*1000):03d}'
    code=f"try, cd('{H.as_posix()}'); run_direct('{case}',{step}); catch ME, disp(getReport(ME)); exit(1); end; exit(0)"
    p=subprocess.run([shutil.which('matlab') or 'matlab','-wait','-singleCompThread','-nosplash','-nodesktop','-minimize','-logfile',str(H/'logs'/(tag+'.log')),'-r',code],creationflags=subprocess.CREATE_NO_WINDOW)
    if p.returncode:raise RuntimeError(tag)
    print(tag,(H/'runs'/tag/'summary.json').read_text(),flush=True)
if __name__=='__main__':
    freeze()
    if '--freeze-only' not in sys.argv:
        for sub in ['logs','runs','results']:(H/sub).mkdir(exist_ok=True)
        with ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(job,[(c,s) for c in ['C1','C2'] for s in [.01,.005]]))
