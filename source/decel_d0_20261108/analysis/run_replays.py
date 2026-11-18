from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,time
H=Path(__file__).resolve().parents[1]
def replay(job):
    m,step=job; tag=f'{m}_{round(step*1000):03d}';folder=H/'runs'/tag
    end=time.monotonic()+1200
    while not (folder/(m+'_batch_summary.json')).exists():
        if time.monotonic()>end:raise TimeoutError(tag)
        time.sleep(2)
    code=f"try, cd('{(H/'analysis').as_posix()}'); replay_d0('{m}',{step}); catch ME, disp(getReport(ME)); exit(1); end; exit(0)"
    cmd=['C:/Program Files/MATLAB/R2018a/bin/matlab.exe','-wait','-singleCompThread','-nosplash','-nodesktop','-minimize','-logfile',str(H/'logs'/f'{tag}_replay.log'),'-r',code]
    subprocess.run(cmd,creationflags=subprocess.CREATE_NO_WINDOW,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    assert (folder/'replay.csv').exists()
    print(tag,flush=True)
if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=10) as ex:list(ex.map(replay,[(m,s) for m in ['O','T','N10','Am','F','D','A'] for s in [.01,.005]]))
