"""Check both QPs on every row and independently enumerate all flagged samples."""
import csv
import itertools
import json
from collections import Counter
from math import degrees

import numpy as np
from scipy.io import loadmat
from check_inputs import HERE, A1, check_inputs

SCALE=np.array([2.,8*np.pi/180,12*np.pi/180])
KEEP=np.r_[0:13,15:18]
EP=[13,14]
TOL=1e-8
KKT=5e-8

def mat(p):return loadmat(p,simplify_cells=True)
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,data):
    with (HERE/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
def write_json(name,obj):
    (HERE/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def assertclose(a,b,tol=1e-11):
    np.testing.assert_allclose(a,b,atol=tol,rtol=3e-12)

def batch_kkt(c,G,u0,u,dual):
    """Recreate centered unit-row KKT with published original-scale duals."""
    raw=-G; row=np.linalg.norm(raw*SCALE,axis=2); nz=row>0
    safe=np.where(nz,row,1.); z=(u-u0)/SCALE
    A=raw*SCALE/safe[:,:,None]
    rhs=(c-np.einsum('nri,ni->nr',raw,u0))/safe
    residual=np.einsum('nri,ni->nr',A,z)-rhs
    mult=dual*row
    azabs=np.einsum('nri,ni->nr',abs(A),abs(z))
    # Constant rows are separately checked in the original physical constraints.
    primal=np.max(np.where(nz,np.maximum(residual,0)/(1+abs(rhs)+azabs),0),axis=1)
    dualnorm=np.max(abs(mult),axis=1); atmu=np.einsum('nri,nr->ni',A,mult)
    zinfty=np.max(abs(z),axis=1)
    st=np.max(abs(z+atmu),axis=1)/(1+zinfty+np.max(abs(atmu),axis=1))
    dn=np.max(np.maximum(-mult,0),axis=1)/(1+dualnorm)
    co=np.max(np.where(nz,abs(mult*residual),0),axis=1)/(1+dualnorm*(1+np.max(np.where(nz,abs(rhs),0),axis=1))+np.sum(z*z,axis=1))
    row_tol=np.where(nz,row*(1e-10+1e-12*(1+abs(rhs)+azabs))+
                     64*np.finfo(float).eps*(1+abs(c)+np.einsum('nri,ni->nr',abs(raw),abs(u))),0)
    return np.maximum.reduce([primal,st,dn,co]), row_tol

def enum_svd(c,G,u0):
    """Separate implementation. Enumerate 0..3 independent normals using SVD."""
    raw=-G*SCALE
    norms=np.linalg.norm(raw,axis=1); nz=norms>0
    assert np.all(c[~nz]>=0)
    A=raw[nz]/norms[nz,None]; b=(c[nz]+G[nz]@u0)/norms[nz]
    best=None; bestval=float('inf'); tried=0
    for size in range(4):
        for support in itertools.combinations(range(len(b)),size):
            if size:
                idx=list(support); active=A[idx]
                U,sv,Vt=np.linalg.svd(active,full_matrices=False)
                if sv[-1]<=1e-12*sv[0]:continue
                coeff=(U.T@b[idx])/sv
                z=Vt.T@coeff; mu=U@(-coeff/sv)
                if min(mu)<-1e-10*(1+np.max(abs(z))):continue
                if np.linalg.norm(z+active.T@mu,np.inf)>1e-8*(1+np.max(abs(z))):continue
            else: z=np.zeros(3)
            tried+=1
            if np.any(A@z-b>1e-10+1e-12*(1+abs(b)+abs(A)@abs(z))):continue
            j=.5*float(z@z)
            if j<bestval:bestval=j;best=u0+SCALE*z
    assert best is not None
    return best,tried

def main():
    assert check_inputs()==json.loads((HERE/'input_hashes.json').read_bytes())
    ev0={(int(r['id']),int(r['matlab_index'])):r for r in read(A1/'event_checks.csv')}
    elig_ids=sorted({k[0] for k in ev0})
    ledger=[];events=[];cases=[];excluded=[];enums=[];maxks=[0.,0.]
    for id in elig_ids:
        tag=f'id_{id:03d}'
        Ymat=mat(A1/'runs'/(tag+'.mat'))['S']; Y=Ymat['Y']
        R=mat(A1/'runs'/(tag+'_rows.mat')); m=mat(HERE/'runs'/(tag+'.mat'));t=m['t']
        X=dict(zip(m['fields'],m['M'].T));status=m['status']
        assert np.array_equal(t,Ymat['t']) and np.array_equal(m['keep']-1,KEEP) and int(m['id'])==id
        eligible=X['base_ok']==1
        valid=np.isfinite(m['u16']).all(axis=1)&np.isfinite(m['u18']).all(axis=1)
        if not valid.all():
            invalid=[(k+1,str(status[k])) for k in np.flatnonzero(~valid)]
            raise ValueError(f'{id}: incomplete solver records; inspect before issuing conclusions: {invalid}')
        n=len(t);assert n==2401
        c,G=R['affineC'],R['affineG']; u0,u,v=m['u0'],m['u18'],m['u16']
        assertclose(Y['u0'],u0);assertclose(Y['uGov'],u)
        kfull,_=batch_kkt(c,G,u0,u,R['dual'])
        kdel,_=batch_kkt(c[:,KEEP],G[:,KEEP,:],u0,v,m['deletion_dual'])
        maxks=[max(maxks[0],max(kfull)),max(maxks[1],max(kdel))]
        assert np.all(kfull[eligible]<=KKT) and np.all(kdel[eligible]<=KKT)
        p18=c+np.einsum('nri,ni->nr',G,u);p16=c+np.einsum('nri,ni->nr',G,v)
        den18=1+abs(c)+np.einsum('nri,ni->nr',abs(G),abs(u))
        den16=1+abs(c)+np.einsum('nri,ni->nr',abs(G),abs(v))
        vn=np.maximum(0,-p16/den16)
        raw=G[:,EP,:];r=np.linalg.norm(raw*SCALE,axis=2);safe=np.where(r==0,1,r)
        z=(v-u0)/SCALE
        b=(c[:,EP]+np.einsum('nri,ni->nr',raw,u0))/safe
        absA=abs(raw*SCALE)/safe[:,:,None]
        tol=np.where(r>0,r*(1e-10+1e-12*(1+abs(b)+np.einsum('nri,ni->nr',absA,abs(z))))+
                     64*np.finfo(float).eps*(1+abs(c[:,EP])+np.einsum('nri,ni->nr',abs(raw),abs(v))),0)
        delta=v-u;d=np.max(abs(delta/SCALE),axis=1)
        j18=.5*np.sum(((u-u0)/SCALE)**2,axis=1);j16=.5*np.sum(z**2,axis=1)
        nested=j16<=j18+TOL*(1+j18)
        row_test=-p16[:,EP]>tol;norm_test=vn[:,EP]>TOL
        confirmed=(d>TOL)&np.any(row_test&norm_test,axis=1)&eligible
        assert np.all(nested[eligible]) and np.all(np.maximum(0,-p18/den18)[eligible]<=TOL) and np.all(vn[eligible][:,KEEP]<=TOL)
        assertclose(p18,Y['phi'],1e-9);assertclose(p16[:,EP],np.stack([X['phi16_lower'],X['phi16_upper']],axis=1),1e-9)
        assertclose(tol,np.stack([X['tol_lower'],X['tol_upper']],axis=1))
        assertclose(vn[:,EP],np.stack([X['norm_lower'],X['norm_upper']],axis=1))
        assertclose(d,X['d']);assertclose(j16,X['objective16']);assertclose(j18,X['objective18'])
        assert np.array_equal(confirmed,X['accepted']==1) and np.all(X['enum_verified'][confirmed]==1)
        for k in np.flatnonzero(X['check_enum']==1):
            a,na=enum_svd(c[k],G[k],u0[k]);b,nb=enum_svd(c[k,KEEP],G[k,KEEP,:],u0[k])
            da=float(max(abs((a-u[k])/SCALE))); db=float(max(abs((b-v[k])/SCALE)))
            assert max(da,db)<=TOL
            ee=c[k,EP]+G[k,EP]@b
            # Confirm the same endpoint exceeds *both* tolerances.
            independent=(max(abs((a-b)/SCALE))>TOL and np.any((-ee>tol[k])&(-ee/(1+abs(c[k,EP])+abs(G[k,EP])@abs(b))>TOL)))
            assert bool(confirmed[k])==bool(independent) and X['enum_verified'][k]==1
            enums.append(dict(id=id,index=k+1,time_s=float(t[k]),matlab_diff18=X['enum_d18'][k],
                              matlab_diff16=X['enum_d16'][k],svd_diff18=da,svd_diff16=db,
                              enum18_kkt=X['enum_kkt18'][k],enum16_kkt=X['enum_kkt16'][k],
                              svd_dual_admissible_candidates18=na,svd_dual_admissible_candidates16=nb))
        seen=[]
        for k in range(n):
            src=ev0.get((id,k+1)); sup=bool(X['support_upper'][k]);assert sup==(src is not None)
            row=dict(id=id,index=k+1,time_s=round(float(t[k]),2),status=str(status[k]),support_upper=int(sup),
                     backup18=int(X['support_backup'][k]),backup16=int(X['retained_backup_support'][k]),
                     d_scaled_inf=d[k],delta_V_mps3=delta[k,0],delta_gamma_degps2=degrees(delta[k,1]),
                     delta_chi_degps2=degrees(delta[k,2]),phi18_lower_Nps=p18[k,13],phi18_upper_Nps=p18[k,14],
                     phi16_lower_Nps=p16[k,13],phi16_upper_Nps=p16[k,14],norm_viol_lower=vn[k,13],norm_viol_upper=vn[k,14],
                     row_tol_lower_Nps=tol[k,0],row_tol_upper_Nps=tol[k,1],objective18=j18[k],objective16=j16[k],
                     kkt18=kfull[k],kkt16=kdel[k],retained_norm_viol=np.max(vn[k,KEEP]),
                     eligible=int(eligible[k]),enumerated=int(X['check_enum'][k]==1),confirmed=int(confirmed[k]))
            ledger.append(row)
            if not eligible[k] or row['status'] not in {'confirmed','support_only','no_detected_change'}:excluded.append(row)
            if sup or row['enumerated']:
                ev=dict(row)
                for name in ['J_included','E_included','attitude_within','backup_VB_over_rho','ref_upper_N','fb_upper_N','Traw_N','Tact_N']:
                    ev[name]=src[name] if src else ''
                events.append(ev);seen.append(ev)
        cases.append(dict(id=id,grid=n,eligible=int(sum(eligible)),upper_support=sum(r['support_upper'] for r in seen),
                         confirmed=int(sum(confirmed)),confirmed_initial=int(confirmed[0]),confirmed_later=int(sum(confirmed[1:])),
                         first_confirmed_s=next((round(float(t[k]),2) for k in np.flatnonzero(confirmed)),None),
                         max_d=float(max(d)),max_kkt18=float(max(kfull)),max_kkt16=float(max(kdel)),
                         max_u0_replay=float(max(X['u0_match'])),max_u18_replay=float(max(X['u18_match'])),max_phi_replay=float(max(X['phi_match']))))
    assert len(ledger)==40817 and len(ev0)==21
    write('all_samples.csv',ledger);write('trial_summary.csv',cases);write('event_comparison.csv',events)
    if enums:write('enumeration_checks.csv',enums)
    with (HERE/'exceptions.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=ledger[0].keys());w.writeheader();w.writerows(excluded)
    hits=[r for r in ledger if r['confirmed']]; next_ids=sorted({r['id'] for r in hits})[:3]
    out=dict(trials=len(cases),total=len(ledger),eligible=sum(r['eligible'] for r in ledger),exceptions=len(excluded),
             status_counts=dict(Counter(r['status'] for r in ledger)),
             exactly_equal=int(sum(r['d_scaled_inf']==0 for r in ledger)),below_change_gate=int(sum(r['d_scaled_inf']<=TOL for r in ledger)),
             confirmed=len(hits),new_outside_a1=sum(not r['support_upper'] for r in hits),enumerated=len(enums),max_independent_kkt=maxks,
             max_solver_comparison=max([max(r['matlab_diff18'],r['matlab_diff16'],r['svd_diff18'],r['svd_diff16']) for r in enums],default=0),
             a3_ids=next_ids,delta_is='u16 minus u18',max_original_replay={key:max(c[key] for c in cases) for key in ('max_u0_replay','max_u18_replay','max_phi_replay')},
             eligible_nonconfirmed_max_d=max((r['d_scaled_inf'] for r in ledger if r['eligible'] and not r['confirmed']),default=0))
    for name,sub in [('initial',[r for r in hits if r['time_s']==0]),('later',[r for r in hits if r['time_s']>0])]:
        objs=[r for r in events if r['confirmed'] and (r['time_s']==0)==(name=='initial')]
        out[name]=dict(points=len(sub),trial_ids=sorted({t['id'] for t in sub}),backup_support=sum(r['backup18'] for r in sub),
                       J_included=sum(int(r['J_included'] or 0) for r in objs),
                       time_range=[min(r['time_s'] for r in sub),max(r['time_s'] for r in sub)] if sub else [],
                       lower_viol=int(sum(r['norm_viol_lower']>TOL and -r['phi16_lower_Nps']>r['row_tol_lower_Nps'] for r in sub)),
                       upper_viol=int(sum(r['norm_viol_upper']>TOL and -r['phi16_upper_Nps']>r['row_tol_upper_Nps'] for r in sub)),
                       d_range=[min(r['d_scaled_inf'] for r in sub),max(r['d_scaled_inf'] for r in sub)] if sub else [],
                       upper_phi_range=[min(r['phi16_upper_Nps'] for r in sub),max(r['phi16_upper_Nps'] for r in sub)] if sub else [])
    write_json('a2_verified.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
