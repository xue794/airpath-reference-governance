"""Disclosed finite follow-up: recentered straight trim and relative turn trim.

The copied kernel is transformed in memory with assertions. Original files are
not changed. For turn trim, delta_nuChi is the fifth Lyapunov coordinate;
physical nuChi=delta_nuChi+nuChiB, uChi=-1.5*delta_nuChi. This preserves A_B.
"""
from pathlib import Path
import types,json,hashlib,sys
HERE=Path(__file__).resolve().parent
source=(HERE/'original_certificate/interval_certificate.py').read_text(encoding='utf-8')
replacements={
 'dv,ga,nv,ng,nc=z;V=Vb+dv': 'dv,ga,nv,ng,dc=z;nc=dc+TURN;V=Vb+dv',
 "uC=-I('1.5')*nc": "uC=-I('1.5')*dc",
 "'phinuCmin':I('1.5')*12*DEG,'phinuCmax':I('1.5')*12*DEG":
 "'phinuCmin':I('1.5')*(12*DEG+TURN),'phinuCmax':I('1.5')*(12*DEG-TURN)"}
for old,new in replacements.items():
    assert source.count(old)==1,old
    source=source.replace(old,new)
source+="\nTURN=I(0)\n"
(HERE/'relative_turn_interval_kernel.py').write_text(source,encoding='utf-8')
mod=types.ModuleType('relative_turn_interval_kernel');exec(compile(source,'relative_turn_interval_kernel.py','exec'),mod.__dict__)
experiments=[('31','0','.002'),('31','0','.003'),('34','.15','.001'),('34','.15','.005')]
results=[]
for center,turn,rho in experiments:
    mod.Vb=mod.I(center);mod.TURN=mod.I(turn)
    z=[mod.I(0) for _ in range(5)];checks,valid=mod.evaluate(z)
    assert valid and all(v.lo>0 for v in checks.values())
    result=mod.verify(rho,maxnodes=10001)
    result.update(center_mps=center,course_rate_center_rad_s=turn,
        maxnodes=10001,all_other_gains_unchanged=True,
        scope='Nominal reference relative-equilibrium ellipsoid; no aircraft certificate or switching theorem.',
        equilibrium_checks={k:v.pair() for k,v in checks.items()},
        backup_row_analytic_lower=(mod.I(rho)*mod.I('.145')).pair(),
        kernel_sha256=hashlib.sha256(source.encode()).hexdigest())
    results.append(result)
    (HERE/'followup_certificates.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps({k:result.get(k) for k in ('center_mps','course_rate_center_rad_s','rho','certified','nodes','seconds','failed')}),flush=True)
