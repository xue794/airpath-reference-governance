"""Compare reconstructed results with archived summaries and manuscript counts."""
from pathlib import Path
import csv,json,math
ROOT=Path(__file__).resolve().parents[1]
def rows(name):
    with (ROOT/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
count=0;maximum=0.
def close(a,b,tol=1e-8):
    global count,maximum
    a,b=float(a),float(b);count+=1
    if math.isnan(a) and math.isnan(b):return
    e=abs(a-b);maximum=max(maximum,e)
    assert e<=tol,(a,b,e,tol)
s=json.loads((ROOT/'generated/summary.json').read_text())
m=s['monte_carlo']
for k,v in {'pairs':500,'both':408,'scalar_only':23,'qp_only':23,'neither':46,'scalar_completed':431,'qp_completed':431,'scalar_lower_rms_common':382}.items():close(m[k],v,0)
assert round(m['median_qp_minus_scalar_rms_m'],3)==3.994
expected={int(r['Index']):r for r in rows('data/original_comparison/tables/W2E_pair_results.csv')}
for r in rows('generated/paired_rms.csv'):
    e=expected[int(r['Index'])]
    close(r['BothCompleted'],e['BothCompleted'],0)
    close(r['dPosRMSFull_m'],e['dPosRMSFull_m'])
expected={(r['Method'],r['Case']):r for r in rows('data/ordinary/case_summary.csv')}
for r in rows('generated/ordinary_metrics.csv'):
    e=expected[r['Method'],r['Case']]
    for k in ['Samples','PositionRMS_grid_m','PositionPeak_m','BackupOutsideSamples','ReferenceThrustMin_N','RawThrustMin_N','RawThrustMax_N','ThrustTargetClippedSamples']:close(r[k],e[k])
    if r['Method']=='Transparent':
        for k,vals in [('reference_below_3N',{'hold':0,'varying':155,'noisy':157}),('reference_negative',{'hold':0,'varying':0,'noisy':17})]:close(r[k],vals[r['Case']],0)
a=s['ablation']; close(a['samples'],40817,0);close(a['trials'],17,0);close(a['confirmed_changes'],21,0)
expected={(r['case'],r['step']):r for r in rows('data/S11/S11_5_direct_center/speed_comparison.csv')}
for r in rows('generated/direct_comparison.csv'):
    e=expected[r['case'],r['step']]
    for k in ['O','Direct','F','O_minus_F','Direct_minus_F','fraction_of_O_F_gain_at_Direct']:close(r[k],e[k])
    assert r['F_extra_gain_below_research_scale']=='True'
assert round(s['timing']['scalar_median_ms'],3)==.778
assert round(s['timing']['qp_median_ms'],3)==1.015
report={'passed':True,'numeric_comparisons':count,'max_absolute_difference':maximum,'scope':'Selected saved-output summaries only; no MATLAB reintegration or interval proof regenerated.'}
(ROOT/'generated/verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('PASS',json.dumps(report))
