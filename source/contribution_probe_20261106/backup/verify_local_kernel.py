"""Validate the three explicit turn-center changes against original and MATLAB."""
from pathlib import Path
import sys,json,importlib.util
from fractions import Fraction
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'original_certificate'))
import interval_certificate as original
spec=importlib.util.spec_from_file_location('turn_kernel',HERE/'relative_turn_interval_kernel.py')
turn=importlib.util.module_from_spec(spec);spec.loader.exec_module(turn)
zero_tests=0
for center in ['31','32','34','38.1']:
    original.Vb=original.I(center);turn.Vb=turn.I(center);turn.TURN=turn.I(0)
    for axis in range(5):
        vals=['0']*5;vals[axis]='.0003'
        a,va=original.evaluate([original.I(v) for v in vals]);b,vb=turn.evaluate([turn.I(v) for v in vals])
        assert va==vb and a.keys()==b.keys()
        assert all((a[k].lo,a[k].hi)==(b[k].lo,b[k].hi) for k in a)
        zero_tests+=1
rows=json.loads((HERE/'turn_matlab_crosscheck.json').read_text(encoding='utf-8'))
keys=['phiVmin','phiVmax','phiGmin','phiGmax','phinuVmin','phinuVmax',
      'phinuGmin','phinuGmax','phinuCmin','phinuCmax','phiAG','phiMuPlus',
      'phiMuMinus','phiNminus','phiNplus','phiTmin','phiTmax']
turn.Vb=turn.I('34');turn.TURN=turn.I('.15');details=[];max_discrepancy=0.
for j,row in enumerate(rows):
    z=[turn.I(str(Fraction.from_float(float(v)))) for v in row['z']]
    bounds,valid=turn.evaluate(z);assert valid
    errors=[]
    for key,value in zip(keys,row['phi']):
        lo,hi=bounds[key].pair();errors.append(max(lo-value,0.,value-hi))
    for key,value in [('hNminus',row['hNminus']),('hNplus',row['hNplus'])]:
        lo,hi=bounds[key].pair();errors.append(max(lo-value,0.,value-hi))
    lo,hi=(bounds['hTmin']+3).pair();errors.append(max(lo-row['Tff'],0.,row['Tff']-hi))
    discrepancy=max(errors);assert discrepancy<1e-8
    max_discrepancy=max(max_discrepancy,discrepancy)
    details.append(dict(point=j,max_MATLAB_distance_outside_interval=discrepancy))
result=dict(zero_turn_original_kernel_bitwise_checks=zero_tests,
    MATLAB_relative_turn_points=len(rows),max_MATLAB_distance_outside_interval=max_discrepancy,
    tolerance=1e-8,passed=True,details=details,
    interpretation='Finite transcription regression only; the separately saved outward interval coverage supplies the all-state nominal reference certificate.')
(HERE/'kernel_regression.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
