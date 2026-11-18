from run_all import run,H
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
if __name__=='__main__':
    tags=json.loads((H/'audit_preflight/restart_jobs.json').read_bytes())
    jobs=[(c,m,int(i)/1000) for c,m,i in (t.split('_') for t in tags)]
    with ThreadPoolExecutor(max_workers=18) as pool:
        for f in as_completed([pool.submit(run,j) for j in jobs]):print(f.result(),flush=True)
