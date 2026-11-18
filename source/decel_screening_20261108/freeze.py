"""Freeze section 2 and its historical inputs before evaluating a new case."""
from pathlib import Path
import hashlib, json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OLD = ROOT / 'endpoint_probe_20261106/full6dof/V31_C60_warm_fine'
AUDIT = ROOT / 'contribution_probe_20261106/backup/V31_backup_saved_audit.mat'
FB = ROOT / 'contribution_probe_20261106/feedback'
CORE = ROOT / 'revision_probe_20261106/source_snapshot/core'
RESULT = HERE / 'results'
RESULT.mkdir(exist_ok=True)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def lock():
    files = [next(ROOT.glob('Aerospace_有限时长*docx')), ROOT/'revision_stage1_20261107/locked_configuration.json',
             AUDIT, ROOT/'revision_wp3_20261107/interval_wp3.py', ROOT/'revision_wp3_20261107/inputs/protocol.json',
             ROOT/'endpoint_probe_20261106/reference/screening_summary.csv',OLD.parent/'ep6_reference.m',
             OLD.parent/'run_ep6_warm.m', FB/'feedback_transfer_audit.py']
    for name in ['T','N10']:
        files.extend([OLD/(name+'.mat'), OLD/(name+'_summary.json'), FB/f'V31_C60_warm_fine_{name}_transfer.csv'])
    files.extend(CORE/f for f in ['rg_constraint_affine.m','rg_N_map.m','rg_alpha_root.m','rg_nominal_aero.m','rg_backup_control.m','rg_backup_state.m'])
    config = dict(stage='Section 2, offline screening only',date='2026-11-08',command=dict(initial=38.1,target_mps=[31,31.5],decel_s=[8,10],start=3,horizon=35),held_out=['C1 31/10 calm','C2 31.5/10 calm','C3 31.5/10 frozen wind'],old=['H31/8/60deg N10','H31/8/60deg T through event'],windows=[[3,16],[14,20],[25,35]],
      center=dict(low=31.25,descent_s=[18,22,26],recovery=[35,55],rho=[.03,.01,.003],sample_dt=.05,sections='At 0,3,5.6,14,14.65,16,16.8,18.7,20,35,40,45,50,55 and at 1/4,1/2,3/4,end of each descent; +/- each Cholesky principal coordinate at 0.97 of radius, 10 points/section.'),
      feedback=dict(direction=[.125,.45,.35],units=['m/s','deg','deg'],uniform='smallest legacy scaled box containing the directional box, scale 38.1',reference_states='All saved T and N10 states and 0.05 s candidate centers in [14,20]. Report rigorous bounds at fixed points, use 8 error vertices as lower bounds on the worst loss.'),
      lp=dict(rows=list(range(1,18)),oldBackupAdmissionTol=1e-9,method='HiGHS max common scaled residual; three free inputs, no invented amplitude bounds; scaled feasibility tol 1e-9',rawAuditTol=1e-8),
      gate='Center regularity, strictly positive operating and state margins plus nonzero tested sections, and plausible feedback bracket. Strict counterexample rejects only the identified candidate. Distinguish 10 N payment from bracket using additional nominal margin. Do not run a regional tube proof or a new trajectory.',
      hashes={str(p.relative_to(ROOT)):digest(p) for p in files})
    p=HERE/'protocol_lock.json'
    if p.exists():
        assert json.loads(p.read_text(encoding='utf8'))==config, 'Frozen input or protocol changed'
    else:
        p.write_text(json.dumps(config,indent=2,ensure_ascii=False),encoding='utf8')

if __name__=='__main__': lock()
