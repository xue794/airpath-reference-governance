"""Recompute selected manuscript summaries from saved outputs, without MATLAB."""
from pathlib import Path
import csv,json,math,statistics
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'generated';OUT.mkdir(exist_ok=True)
def rows(name):
    with (ROOT/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def number(x):return float(x)
def truth(x):return str(x).lower() in ['1','true']
def save(name,data):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def mean_abs_speed(d):
    w=[r for r in d if 3-1e-9<=float(r['time'])<=16+1e-9]
    assert abs(float(w[0]['time'])-3)<1e-9 and abs(float(w[-1]['time'])-16)<1e-9
    y=[abs(float(r['vground'])-float(r['vrequest'])) for r in w]
    return math.fsum((float(b['time'])-float(a['time']))*(u+v)/2 for a,b,u,v in zip(w,w[1:],y,y[1:]))/13

# Derive pair membership and differences from per-method records, not pair summaries.
methods=rows('data/original_comparison/tables/W2E_method_results.csv')
by={m:{int(r['Index']):r for r in methods if r['Method']==m} for m in ['Proposed','MatchedQP']}
assert set(by['Proposed'])==set(by['MatchedQP'])==set(range(1,501))
pairs=[];counts=Counter()
for i in range(1,501):
    a,b=by['Proposed'][i],by['MatchedQP'][i];ac,bc=truth(a['Completed']),truth(b['Completed'])
    counts['both' if ac and bc else 'scalar_only' if ac else 'qp_only' if bc else 'neither']+=1
    pairs.append({'Index':i,'BothCompleted':int(ac and bc),'dPosRMSFull_m':float(b['PosRMSFull_m'])-float(a['PosRMSFull_m']) if ac and bc else 'NaN'})
common=[r for r in pairs if r['BothCompleted']]
mc={'pairs':500,**counts,'scalar_completed':sum(truth(x['Completed']) for x in by['Proposed'].values()),
    'qp_completed':sum(truth(x['Completed']) for x in by['MatchedQP'].values()),
    'scalar_lower_rms_common':sum(r['dPosRMSFull_m']>0 for r in common),
    'median_qp_minus_scalar_rms_m':statistics.median(r['dPosRMSFull_m'] for r in common),
    'scalar_median_rms_m':statistics.median(float(by['Proposed'][r['Index']]['PosRMSFull_m']) for r in common),
    'qp_median_rms_m':statistics.median(float(by['MatchedQP'][r['Index']]['PosRMSFull_m']) for r in common)}
save('paired_rms.csv',pairs)

ordinary=[]
for method in ['Transparent','Proposed','MatchedQP']:
    for case in ['hold','varying','noisy']:
        d=rows(f'data/ordinary/{method}_{case}_layered_timeseries.csv')
        v=lambda k:[float(r[k]) for r in d]
        ordinary.append({'Method':method,'Case':case,'Samples':len(d),'PositionRMS_grid_m':math.sqrt(math.fsum(x*x for x in v('position_error_m'))/len(d)),
            'PositionPeak_m':max(v('position_error_m')),'reference_below_3N':sum(x<3 for x in v('reference_T_N')),
            'reference_negative':sum(x<0 for x in v('reference_T_N')),'BackupOutsideSamples':sum(x>1+1e-8 for x in v('VB_over_rho')),
            'ReferenceThrustMin_N':min(v('reference_T_N')),'RawThrustMin_N':min(v('raw_T_N')),'RawThrustMax_N':max(v('raw_T_N')),
            'ThrustTargetClippedSamples':sum(truth(r['thrust_target_clipped']) for r in d)})
save('ordinary_metrics.csv',ordinary)
ab=rows('data/endpoint_extension/a2/all_samples.csv')
ablation={'samples':len(ab),'trials':len({r['id'] for r in ab}),'confirmed_changes':sum(truth(r['confirmed']) for r in ab)}
timing=rows('data/value_checks/curve/timing_paired_blocks.csv')
s=statistics.median(float(r['whole_scalar_ms']) for r in timing);q=statistics.median(float(r['whole_qp_ms']) for r in timing)
cost={'blocks':len(timing),'scalar_median_ms':s,'qp_median_ms':q,'reduction_ratio_of_medians_percent':100*(q-s)/q}
direct=[]
for case in ['C1','C2']:
    d={arm:rows(f'data/direct_traces/{case}_{arm}_005.csv') for arm in ['O','Direct','F']}
    m={arm:mean_abs_speed(trace) for arm,trace in d.items()}
    gap=lambda arm:min(float(r['fbHi']) for r in d[arm] if 14-1e-9<=float(r['time'])<=20+1e-9)
    rms=lambda arm:math.sqrt(math.fsum(float(r['positionError'])**2 for r in d[arm])/len(d[arm]))
    direct.append({'case':case,'step':'005',**m,'O_minus_F':m['O']-m['F'],'Direct_minus_F':m['Direct']-m['F'],
        'fraction_of_O_F_gain_at_Direct':(m['O']-m['Direct'])/(m['O']-m['F']),
        'position_rms_Direct_minus_F_m':rms('Direct')-rms('F'),
        'minimum_upper_gap_Direct_minus_F_N':gap('Direct')-gap('F'),
        'F_extra_gain_below_research_scale':m['Direct']-m['F']<.606})
save('direct_comparison.csv',direct)
(OUT/'summary.json').write_text(json.dumps({'scope':'Saved-output reconstruction; no simulation or interval-tree rerun','monte_carlo':mc,'ablation':ablation,'timing':cost,'direct':direct},indent=2),encoding='utf-8')
print('Reconstructed paired, ordinary, ablation, timing and direct-center summaries in generated/.')
