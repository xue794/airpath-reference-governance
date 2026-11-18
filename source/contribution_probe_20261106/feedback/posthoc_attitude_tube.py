"""Supplement only: conditional raw-thrust tube at theta <= 2.5 degrees.

2.5 degrees was selected AFTER seeing the largest attitude error (2.325 deg)
in the two archived N10 runs. It is NOT a proven closed-loop bound.
The 0.70-degree columns and original complete interval audit are retained.
"""
import csv,json
from concurrent.futures import ProcessPoolExecutor
from feedback_transfer_audit import ROOT,SOURCE,loadmat,prepare,radius

def run(case):
    L=loadmat(SOURCE/case/'N10.mat',simplify_cells=True);P=prepare(L['C'])
    with (ROOT/f'{case}_N10_transfer.csv').open(newline='',encoding='utf-8') as f:old=list(csv.DictReader(f))
    rows=[]
    for k,(x,r) in enumerate(zip(L['Xs'],old)):
        rad=radius(x[19:25],P,'physical',theta_deg=2.5)
        err=float(r['error_box_mps']);th=float(r['attitude_error_deg'])
        rows.append(dict(time_s=float(r['time_s']),error_box_mps=err,attitude_error_deg=th,
            radius_four_physical_theta2p5_mps=rad,physical_theta2p5_inclusion=bool(rad>0 and err<=rad and th<=2.5)))
        if k%100==0:print(case,k,len(old),flush=True)
    with (ROOT/f'{case}_N10_theta2p5.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    return dict(case=case,samples=len(rows),theta_deg=2.5,selection='post hoc rounded-up diagnostic assumption; not a prior theorem',
        minimum_radius_mps=min(r['radius_four_physical_theta2p5_mps'] for r in rows),
        inclusion_samples=sum(r['physical_theta2p5_inclusion'] for r in rows),
        attitude_condition_samples=sum(r['attitude_error_deg']<=2.5 for r in rows))

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=2) as pool:result=list(pool.map(run,['V31_C60_warm_fine','V32_C60_warm_fine']))
    (ROOT/'posthoc_attitude_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('POSTHOC_ATTITUDE_AUDIT_COMPLETE',flush=True)
