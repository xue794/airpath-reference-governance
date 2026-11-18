"""Check the representative C1 Direct MATLAB rerun against the archived trace."""
from pathlib import Path
import argparse,csv,json,math

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--run-folder',type=Path,default=ROOT/'generated/matlab_C1_005')
args=parser.parse_args()
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
folder=args.run_folder/'source/decel_center_20261108/runs/C1_005'
actual=rows(folder/'trace.csv');expected=rows(ROOT/'data/direct_traces/C1_Direct_005.csv')
assert len(actual)==len(expected)==701
assert list(actual[0])==list(expected[0])
summary=json.loads((folder/'summary.json').read_text())
assert summary['status']=='completed' and summary['end_time']==35 and not summary['events_s'] and not summary['event_indices']
max_by_field={k:0. for k in actual[0]};comparisons=0
for a,b in zip(actual,expected):
    for key in a:
        x,y=float(a[key]),float(b[key]);comparisons+=1
        if math.isnan(x) and math.isnan(y):continue
        assert math.isfinite(x) and math.isfinite(y),(key,x,y)
        delta=abs(x-y);max_by_field[key]=max(max_by_field[key],delta)
        assert delta<=1e-8,(key,a['time'],x,y,delta)
def metrics(d):
    q=[r for r in d if 3-1e-9<=float(r['time'])<=16+1e-9]
    errors=[abs(float(r['vground'])-float(r['vrequest'])) for r in q]
    mae=math.fsum((float(b['time'])-float(a['time']))*(x+y)/2 for a,b,x,y in zip(q,q[1:],errors,errors[1:]))/13
    rms=math.sqrt(math.fsum(float(r['positionError'])**2 for r in d)/len(d))
    upper=min(float(r['fbHi']) for r in d if 14-1e-9<=float(r['time'])<=20+1e-9)
    return {'ground_speed_MAE_3_16_mps':mae,'position_RMS35_m':rms,'min_upper_feedback_gap_14_20_N':upper}
a,b=metrics(actual),metrics(expected)
for key in a:assert abs(a[key]-b[key])<=1e-8,(key,a[key],b[key])
run=json.loads((args.run_folder/'release_run.json').read_text())
report={'passed':True,'case':'C1','arm':'Direct','max_step_s':.005,'horizon_s':35,'samples':701,
    'trace_values_compared':comparisons,'trace_absolute_tolerance':1e-8,
    'max_trace_absolute_difference':max(max_by_field.values()),'max_difference_by_field':max_by_field,
    'reintegrated_metrics':a,'archived_metrics':b,
    'integration_wall_seconds':summary['wall_seconds'],'matlab_version':run['matlab_version'],
    'scope':'One original case reintegrated with unchanged source and parameters in an isolated directory; no all-study validation claim.'}
(args.run_folder/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('PASS',json.dumps(report,indent=2))
