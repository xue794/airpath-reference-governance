"""Native-kernel regression, terminal handoff, and direct uniform-box witnesses."""
from model import *
from kernel import guide,nominal,domain,Box,box
from scipy.io import loadmat
sys.path.insert(0,str(OLD))
from screen_feedback import ForceCheck,guide as exact_guide,BOXES,pt
def run():
    lock()
    assert R['alphaLo']==0 and R['Tmin']==0 and R['Tmax']==60 and R['kT']==R['kN']==1 and R['DeltaN']==10 and R['DeltaT']==3
    assert R['kB']>=.12 and R['gammaB']==0 and R['kChiD']==1.5
    assert np.array_equal(R['AB'],[[0,0,1,0,0],[0,0,0,1,0],[-R['kVp'],0,-R['kVd'],0,0],[0,-R['kGammaP'],0,-R['kGammaD'],0],[0,0,0,0,-R['kChiD']]])
    s=loadmat(OUT/'native_check.mat',simplify_cells=True)
    t=s['T'];x=s['X'];mapping=['Vmin','Vmax','Gmin','Gmax','nuVmin','nuVmax','nuGmin','nuGmax','nuCmin','nuCmax','AG','MuPlus','MuMinus','Nminus','Nplus','Tmin','Tmax']
    phases=np.select([t<=3,t<=25,t<=35],[0,1,2],default=3);maxerr=0.
    for ph in range(4):
        ix=phases==ph;tt=box(t[ix]);v,n,u=guide(tt,22,ph)
        xx=x[ix].T
        z=[box(xx[0])-v,box(xx[1]),box(xx[3])-n,box(xx[4]),box(xx[5])]
        rows=nominal(tt,z,22,ph)
        for j,k in enumerate(mapping):
            b=rows['phi'+k];obs=s['phi'][ix,j]
            maxerr=max(maxerr,float(max(np.max(b.l-obs),np.max(obs-b.h))))
    assert maxerr<1e-9 and s['n']==1271
    calc=ForceCheck(data['C'])
    uniform=calc.examine(exact_guide(14,22),'uniform')
    assert uniform['upper_loss_sup_lb_N']>10
    J=json.loads((OUT/'raw_thrust_0.01.json').read_bytes())
    assert next(v for v in J if v['box']=='uniform')['counterexample']
    groups={}
    for name,parts in [('nominal',['nominal','nominal_focused']),('directional',['directional','directional_focused','joint']),('uniform',['uniform','uniform_focused'])]:
        a=[json.loads((OUT/f'22_0.01_{k}.json').read_bytes()) for k in parts]
        groups[name]=dict(nodes=sum(d['nodes'] for d in a),cpu_s=sum(d['cpu_s'] for d in a))
        assert groups[name]['nodes']<=1000000 and groups[name]['cpu_s']<14400
    checks={}
    for k in ['nominal_focused','directional_focused','uniform_focused','joint']:
        a=json.loads((OUT/f'22_0.01_{k}_review.json').read_bytes())
        checks[k]={row:min(i['worsts'][row]['lower'] for i in a.values()) for row in a['original']['worsts']}
    result=dict(candidate=dict(D=22,rho=.01,center_low=31.25),passed=True,full_reference_window=[0,55],
                feedback_window=[14,20],native_points=int(s['n']),native_max_discrepancies=s['mx'].tolist(),
                interval_to_native_max_difference=maxerr,geometry=geometry(.01),budget=groups,
                whole_tube_bounds=checks,uniform_counterexample_at_center_t14=uniform,
                controller_feedback_witness=J,confirmation_flights=0)
    save(OUT/'decision.json',result)
if __name__=='__main__':run()
