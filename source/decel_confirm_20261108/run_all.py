"""Run the frozen 30 trajectories, then localize chart exits and compare no-tube inputs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib, json, subprocess, shutil
H=Path(__file__).resolve().parent
MATLAB=shutil.which('matlab') or r'C:\Program Files\MATLAB\R2018a\bin\matlab.exe'
def invoke(name, commands):
    code=f"try, cd('{H.as_posix()}'); addpath('analysis'); {commands}; catch ME, disp(getReport(ME)); exit(1); end; exit(0)"
    log=H/'logs'/(name+'.log')
    p=subprocess.run([MATLAB,'-wait','-singleCompThread','-nosplash','-nodesktop','-minimize','-logfile',str(log),'-r',code],creationflags=subprocess.CREATE_NO_WINDOW,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'{name}: MATLAB exit {p.returncode}; {log}')
def run(job):
    sc,m,s=job; tag=f'{sc}_{m}_{round(s*1000):03}'; d=H/'runs'/tag
    args=f"'{sc}','{m}',{s}"
    commands=[]
    if not (d/(m+'_batch_summary.json')).exists(): commands.append('run_confirm('+args+')')
    if not (d/'replay.csv').exists(): commands.append('replay_confirm('+args+')')
    if commands:invoke(tag,'; '.join(commands))
    assert (d/'replay.csv').is_file(),tag
    return tag
if __name__=='__main__':
    cfg=json.loads((H/'protocol.json').read_bytes())
    for path,expected in cfg['inputs'].items(): assert hashlib.sha256((H.parent/path).read_bytes()).hexdigest()==expected,path
    jobs=[(case,m,s) for s in cfg['steps'] for case,arms in cfg['arms'].items() for m in arms]
    failed=[]
    with ThreadPoolExecutor(max_workers=24) as pool:
        futures={pool.submit(run,j):j for j in jobs}
        for f in as_completed(futures):
            try: print('saved',f.result(),flush=True)
            except Exception as e:failed.append(str(e));print('ERROR',e,flush=True)
    if failed: raise RuntimeError(failed)
    invoke('event_free_checks','confirm_chart_events; confirm_free')
    print('ALL CONFIRMATION RUNS AND EVENT CHECKS FINISHED',flush=True)
