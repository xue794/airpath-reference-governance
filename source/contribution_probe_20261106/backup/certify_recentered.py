"""Small disclosed search: exact interval kernel unchanged, only trim Vb changed.

This certifies nominal reference ellipsoids, not six-DOF aircraft trajectories.
An unresolved interval cover is not a counterexample or proof of infeasibility.
"""
from pathlib import Path
import sys, json, hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'original_certificate'))
import interval_certificate as ic
from selftest_interval import main as selftest

selftest()
kernel=HERE/'original_certificate/interval_certificate.py'
experiments=[('31','.001'),('31','.005'),('31','.01'),
             ('32','.01'),('32','.03'),('34','.03'),('34','.05'),('38.1','.03')]
results=[]
for center,rho in experiments:
    ic.Vb=ic.I(center)
    result=ic.verify(rho,maxnodes=10001)
    result.update(center_mps=center, unchanged_feedback_gains=True,
                  kernel_sha256=hashlib.sha256(kernel.read_bytes()).hexdigest(),
                  nominal_reference_only=True, search_maxnodes=10001)
    # Translation leaves A_B,Q_B,P_B and their exact Lyapunov identity intact.
    # Independently verify Q_B-.145 P_B positive definite with the same kernel.
    a,b,c=ic.PV;d,e,f=ic.PG;mu=ic.I('.145')
    checks={'speed_diagonal':ic.q1-mu*a,
            'speed_determinant':(ic.q1-mu*a)*(ic.qv-mu*c)-(mu*b).sq(),
            'climb_diagonal':ic.qg-mu*d,
            'climb_determinant':(ic.qg-mu*d)*(ic.qg-mu*f)-(mu*e).sq(),
            'course_diagonal':ic.qc-mu*ic.PC}
    assert all(z.lo>0 for z in checks.values())
    result['backup_row_analytic']={'mu_lower':'.145','all_PD_tests_pass':True,
        'tests':{k:v.pair() for k,v in checks.items()},
        'strict_lower_bound':(ic.I(rho)*mu).pair(),
        'formula':'rho * min(kB, .145), kB=.4'}
    results.append(result)
    (HERE/'recentered_certificates.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps({k:result.get(k) for k in ('center_mps','rho','certified','nodes','seconds','failed')}),flush=True)

