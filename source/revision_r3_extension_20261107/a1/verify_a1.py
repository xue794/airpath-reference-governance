"""Recompute the A1 reproduction gates and instantaneous tests from the saved data."""
from pathlib import Path
import csv
from collections import Counter
import hashlib
import json
import platform
import sys

import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from scipy import __version__ as scipy_version

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
WP3 = PROJECT / "revision_wp3_20261107"
sys.path.insert(0, str(WP3))
from interval_wp3 import Evaluator, CAP, THETA, SCALE

IDS = [22, 32, 33, 80, 86, 91, 108, 113, 164, 184, 236, 322, 374, 386, 397, 472, 483]
TWO = {22, 86, 108, 374}
CORE = PROJECT / "FWUAV/T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel"
ARCHIVE = CORE.with_name(CORE.name + "_W2E_PairedMC500_rho010_v1_0")
PKG = CORE.parent / "P1_Aerospace_W2E_R2018a_v1_0"

def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def dump(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

def write_csv(name, rows):
    with (HERE / name).open("w",encoding="utf-8-sig",newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def numeric_record(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(numeric_record(a[k],b[k]) for k in a)
    return np.array_equal(a,b)

def maxabs(x): return float(np.max(np.abs(x)))

def verify_sources():
    pkg_hashes=json.loads((PKG/'RELEASE_MANIFEST.json').read_bytes())['files']
    path_list=[PKG/p for p in pkg_hashes]
    assert all(sha(PKG/p)==h for p,h in pkg_hashes.items())
    manifest=read_csv(HERE/'identity/W2E_source_identity.csv')
    assert len(manifest)==66 and all(r['SourceTextMatch']=='1' and sha(CORE/r['File'])==r['ExpectedRawSHA256'] for r in manifest)
    path_list += [CORE/r['File'] for r in manifest]
    path_list += [ARCHIVE/p for p in ['W2E_constraint_role_ledger.csv','W2E_method_results.csv','W2E_CHECKPOINT.mat','W2E_initial_admissibility.csv']]
    path_list += [PKG/'RELEASE_MANIFEST.json', WP3/'interval_wp3.py',WP3/'inputs/protocol.json', PROJECT/'Aerospace_R3端点约束与反馈预留验证方案_v1_0.docx']
    identities=read_csv(HERE/'identity/W2E_sample_identity.csv')[0]
    assert identities['SHA256']==pkg_hashes['data/Rev41Results_FROZEN_SAMPLE_SOURCE.mat'] and identities['ExactRegenerationMatch']=='1'
    assert identities['MismatchCount']=='0' and identities['N']=='500' and identities['Seed']=='20260819'
    dump(HERE/'identity'/'input_hashes.json',{p.relative_to(PROJECT).as_posix():sha(p) for p in path_list})
    return len(pkg_hashes)

def event_checks(id,k,S,rows,V,ga,chi,kappa,Ract,assessor):
    Y=S['Y']; x=S['x']; t=S['t']; xf=x[k,19:25]; C=S['C']; R=C['R']
    dv=V[k]-xf[0]; dg=ga[k]-xf[1]; dc=np.arctan2(np.sin(chi[k]-xf[2]),np.cos(chi[k]-xf[2]))
    e=np.array([dv,dg,dc]); obs=float(max(abs(dv),SCALE*abs(dg),SCALE*abs(dc)))
    theta=float(np.rad2deg(Y['thetaF'][k]+Y['thetaR'][k]))
    assessor.prepare(xf,rows['alphaRef'][k])
    re,ok,_,_,_=assessor.radius('E',CAP)
    rg,_,_,_,_=assessor.radius('G',re)
    rj,okj,_,_,_=assessor.radius('J',rg)
    ref=assessor.numeric(xf,np.zeros(3))[0]
    fb,X,ag,ac,_=assessor.numeric(xf,e)
    af=rows['alphaFeedback'][k]
    cd=np.polynomial.polynomial.polyval(af,assessor.CD)
    raw=(0.5*R['rho']*R['P']['S']*V[k]**2*cd+X)/kappa[k]
    root=rows['alphaRef'][k]; ref_x=assessor.numeric(xf,np.zeros(3))[1]
    ref_drag=.5*R['rho']*R['P']['S']*xf[0]**2*np.polynomial.polynomial.polyval(root,assessor.CD)
    Tff=(ref_drag+ref_x)/np.cos(root)
    assert abs(Tff-(Y['stateMargins'][k,15]+R['Tmin']+R['DeltaT']))<1e-7
    target=float(np.clip(raw,C['P7']['Tmin'],C['P7']['Tmax']))
    np.testing.assert_allclose(ref,rows['gapRef'][k],rtol=0,atol=1e-7)
    np.testing.assert_allclose(fb,rows['gapFeedback'][k],rtol=0,atol=1e-7)
    assert abs(raw-Y['Tcmd'][k])<1e-7
    mu=np.arctan2(ac,ag)
    commanded=(Rotation.from_euler('Z',chi[k])*Rotation.from_euler('Y',ga[k])*
               Rotation.from_euler('X',mu)*Rotation.from_euler('Y',af))
    filtered=Rotation.from_quat(np.r_[x[k,26:29],x[k,25]])
    theta_check=float(np.rad2deg((commanded.inv()*filtered).magnitude()+(filtered.inv()*Ract[k]).magnitude()))
    assert abs(theta-theta_check)<1e-5 and (theta<=THETA)==(theta_check<=THETA)
    primal=float(Y['phi'][k,14])
    return dict(id=id,matlab_index=k+1,time_s=float(t[k]),upper_dual=float(rows['dual'][k,14]),
                upper_dual_scaled=float(rows['dualScaled'][k,14]),upper_phi_Nps=primal,
                backup_support=int(Y['qpSupport'][k,17]),backup_VB_over_rho=float(rows['VBoverRho'][k]),
                eV_mps=dv,eGamma_deg=np.rad2deg(dg),eChi_deg=np.rad2deg(dc),error_box_mps=obs,
                attitude_sum_deg=theta,attitude_rebuilt_deg=theta_check,attitude_within=int(theta<=THETA),
                width_E_mps=re,width_G_mps=rg,width_J_mps=rj,E_included=int(ok and obs<=re),
                G_included=int(rg>0 and obs<=rg),J_air_included=int(okj and obs<=rj),
                J_included=int(okj and obs<=rj and theta<=THETA),
                ref_lower_N=float(ref[0]),ref_upper_N=float(ref[1]),fb_lower_N=float(fb[0]),fb_upper_N=float(fb[1]),
                alpha_ref_deg=float(np.rad2deg(rows['alphaRef'][k])),alpha_fb_deg=float(np.rad2deg(af)),
                Tff_N=float(Tff),Traw_N=float(raw),Ttarget_clipped_N=target,Tact_N=float(Y['Tact'][k]),
                raw_thrust_inside=int(0<raw<60))

def audit_trial(id, old, old_rows, base, sample):
    tag=f'id_{id:03d}'
    chk=json.loads((HERE/'checks'/f'{tag}.json').read_text())
    assert chk['Pass'], f'{tag} did not pass MATLAB reproduction: {chk}'
    load=lambda p: loadmat(p,simplify_cells=True)
    case=load(HERE/'runs'/f'{tag}.mat'); rows=load(HERE/'runs'/f'{tag}_rows.mat')
    S=case['S'];Y=S['Y'];C=S['C'];x=S['x'];t=S['t']
    assert numeric_record(sample,case['sample']) and numeric_record(base['R'],C['R']) and numeric_record(base['A'],C['A'])
    assert case['E']['CorePass']==case['old']['CorePass']==int(old['CorePass'])==1
    np.testing.assert_allclose(t, np.arange(2401)/20, rtol=0, atol=2e-14)
    assert x.shape==(2401,36) and np.isfinite(x).all()
    rd=np.column_stack((C['VgVarying']*t,50*np.sin(.01*t),-C['R']['h_ref']-6*np.sin(.008*t)))
    assert maxabs(rd-Y['rd'])<1e-10
    dist2=np.sum((x[:,:3]-rd)**2,axis=1)
    rms_full=np.sqrt(np.mean(dist2));rms_late=np.sqrt(np.mean(dist2[t>=60]))
    for result,attr in [(rms_full,'PosRMSFull_m'),(rms_late,'PosRMSLate_m')]:
        assert abs(result-case['E'][attr])<1e-10 and abs(result-float(old[attr]))<=max(1e-6,1e-8*abs(float(old[attr])))
    wind=sample['windMean']+sample['windAmp']*np.sin(t[:,None]*sample['windOmega']+sample['windPhase'])
    rot=Rotation.from_quat(x[:,[7,8,9,6]])
    air=rot.apply(x[:,3:6])-wind; body=rot.inv().apply(air);V=np.linalg.norm(air,axis=1)
    ga=np.arctan2(-air[:,2],np.hypot(air[:,0],air[:,1]));chi=np.arctan2(air[:,1],air[:,0]);kappa=body[:,0]/V
    for v,g in [(V,'V'),(ga,'gamma'),(chi,'chi'),(kappa,'kappa'),(wind,'w')]:
        assert maxabs(v-Y[g])<1e-10
    support=np.asarray(Y['qpSupport'],bool); binding=np.asarray(Y['binding'],bool)
    dual=rows['dualScaled']; calc_support=dual>1e-9*(1+np.max(abs(dual),axis=1)[:,None])
    assert np.array_equal(calc_support,support)
    c=rows['affineC'];G=rows['affineG'];u=Y['uGov']
    phi=c+np.einsum('nri,ni->nr',G,u);den=1+abs(c)+np.einsum('nri,ni->nr',abs(G),abs(u))
    assert maxabs(phi-Y['phi'])<1e-7
    assert np.array_equal(np.abs(phi)<=1e-8*den,binding)
    kkt=max(float(np.max(Y[f])) for f in ('qpPrimal','qpStationarity','qpDualViolation','qpComplementarity'))
    normalized=max(0.,float(np.max(-phi/den)))
    xi=Y['xf'][:,[0,1,3,4,5]].copy();xi[:,0]-=C['R']['Vb'];xi[:,1]-=C['R']['gammaB']
    vb=np.einsum('ni,ij,nj->n',xi,C['R']['PB'],xi)/C['R']['rhoB']
    assert maxabs(vb-rows['VBoverRho'])<1e-10 and np.max(vb)<=1+1e-8 and kkt<=5e-8 and normalized<=1e-8
    initial=next(r for r in read_csv(ARCHIVE/'W2E_initial_admissibility.csv') if int(r['Index'])==id)
    assert initial['AllPass']=='1' and initial['QPInitialSolvePass']=='1'
    assert abs(float(initial['RawCommandSpeed_mps'])-np.linalg.norm(Y['vc'][0]))<1e-10
    assert abs(float(initial['InitialVB_over_rhoB'])-vb[0])<1e-12
    assert abs(float(initial['MinBackupSigma'])-np.min(c[0]+G[0]@Y['uB'][0]))<1e-10
    assert abs(float(initial['QPMinPhi'])-min(phi[0]))<1e-10
    roles=[]
    for j,row in enumerate(old_rows):
        assert row['Row']==rows['names'][j] and int(row['RowIndex'])==j+1
        rc=int(sum(support[:,j]));bc=int(sum(binding[:,j]))
        assert rc==int(row['SelectionRoleCount']) and bc==int(row['BindingCount'])
        roles.append(dict(id=id,row=j+1,name=row['Row'],support_archive=int(row['SelectionRoleCount']),support_rerun=rc,
                          binding_archive=int(row['BindingCount']),binding_rerun=bc))
    idx=np.flatnonzero(support[:,14]);assert len(idx)==1+(id in TWO)
    assessor=Evaluator(C)
    events=[event_checks(id,int(k),S,rows,V,ga,chi,kappa,rot,assessor) for k in idx]
    result=dict(id=id,samples=len(t),full_old_m=chk['OldFull'],full_new_m=chk['RMSFull'],full_abs_diff_m=abs(chk['RMSFull']-chk['OldFull']),
                late_old_m=chk['OldLate'],late_new_m=chk['RMSLate'],late_abs_diff_m=abs(chk['RMSLate']-chk['OldLate']),
                upper_support=len(idx),support_times_s=' / '.join(f'{t[k]:.2f}' for k in idx),
                attitude_peak_deg=float(np.rad2deg(np.max(Y['thetaF']+Y['thetaR']))),
                kkt_max=chk['KKT'],norm_row_violation_max=chk['RowViolation'],VB_over_rho_max=chk['MaxVBoverRho'],
                full_grid_python_rms_m=float(rms_full),late_grid_python_rms_m=float(rms_late),python_phi_abs_diff_max=maxabs(phi-Y['phi']),
                integration_s=chk['Integration_s'],core_pass=1,all_18_support_binding_match=1)
    return result,roles,events

def main():
    npkg=verify_sources()
    expected=loadmat(HERE/'identity/frozen.mat',simplify_cells=True)
    original={int(r['Index']):r for r in read_csv(ARCHIVE/'W2E_method_results.csv') if r['Method']=='MatchedQP'}
    ledger=read_csv(ARCHIVE/'W2E_constraint_role_ledger.csv')
    selected=sorted(int(r['Index']) for r in ledger if r['Method']=='MatchedQP' and r['Row']=='N_plus' and int(r['SelectionRoleCount'])>0)
    assert selected==IDS and len(IDS)==len(expected['Samples'])
    allruns=[];allrows=[];events=[]
    for id,sample in zip(IDS,expected['Samples']):
        old_rows=sorted((r for r in ledger if int(r['Index'])==id and r['Method']=='MatchedQP'),key=lambda r:int(r['RowIndex']))
        assert len(old_rows)==18
        result,roles,ev=audit_trial(id,original[id],old_rows,expected['Cbase'],sample)
        allruns.append(result);allrows+=roles;events+=ev
    assert len(events)==21 and len(allrows)==17*18
    write_csv('trial_reproduction.csv',allruns);write_csv('row_counts.csv',allrows);write_csv('event_checks.csv',events)
    summary=dict(phase='A1',accepted=len(allruns),grid_points=sum(r['samples'] for r in allruns),
        upper_support_samples=len(events),role_rows=len(allrows),source_files=66,package_files=npkg,
        rms_full_max_abs_diff_m=max(r['full_abs_diff_m'] for r in allruns),
        rms_late_max_abs_diff_m=max(r['late_abs_diff_m'] for r in allruns),
        min_event_time_s=min(e['time_s'] for e in events),max_event_time_s=max(e['time_s'] for e in events),
        event_time_counts=dict(sorted(Counter(round(e['time_s'],2) for e in events).items())),
        coactive_backup=sum(e['backup_support'] for e in events),attitude_inside=sum(e['attitude_within'] for e in events),
        E_inside=sum(e['E_included'] for e in events),G_inside=sum(e['G_included'] for e in events),
        J_air_inside=sum(e['J_air_included'] for e in events),J_inside=sum(e['J_included'] for e in events),
        fb_root_bracket_observed=sum(min(e['fb_lower_N'],e['fb_upper_N'])>0 for e in events),
        raw_thrust_inside=sum(e['raw_thrust_inside'] for e in events),
        event_error_range_mps=[min(e['error_box_mps'] for e in events),max(e['error_box_mps'] for e in events)],
        event_attitude_range_deg=[min(e['attitude_sum_deg'] for e in events),max(e['attitude_sum_deg'] for e in events)],
        event_fb_upper_range_N=[min(e['fb_upper_N'] for e in events),max(e['fb_upper_N'] for e in events)],
        event_ref_upper_range_N=[min(e['ref_upper_N'] for e in events),max(e['ref_upper_N'] for e in events)],
        max_abs_upper_phi_Nps=max(abs(e['upper_phi_Nps']) for e in events),
        event_raw_thrust_range_N=[min(e['Traw_N'] for e in events),max(e['Traw_N'] for e in events)],
        python=platform.python_version(),numpy=np.__version__,scipy=scipy_version)
    for name, group in [('initial_events',[e for e in events if e['time_s']==0]),
                        ('later_events',[e for e in events if e['time_s']>0])]:
        summary[name]=dict(samples=len(group),trial_ids=sorted({e['id'] for e in group}),
                           time_s=[float(min(e['time_s'] for e in group)),float(max(e['time_s'] for e in group))] if group else None,
                           backup_support=sum(e['backup_support'] for e in group),E_inside=sum(e['E_included'] for e in group),
                           G_inside=sum(e['G_included'] for e in group),J_inside=sum(e['J_included'] for e in group),
                           attitude_inside=sum(e['attitude_within'] for e in group))
    summary['J_exclusions']=dict(air_only=sum(not e['J_air_included'] and e['attitude_within'] for e in events),
                                 attitude_only=sum(e['J_air_included'] and not e['attitude_within'] for e in events),
                                 both=sum(not e['J_air_included'] and not e['attitude_within'] for e in events))
    dump(HERE/'a1_verified.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
