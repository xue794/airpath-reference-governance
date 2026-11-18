from feedback_transfer_audit import *

def attribution(xf,e,P):
    h0=calc(xf,np.zeros(3),P)[0][1]
    def loss(z):return h0-calc(xf,z,P)[0][1]
    out=np.zeros(3)
    for order in itertools.permutations(range(3)):
        z=np.zeros(3);previous=loss(z)
        for i in order:
            z[i]=e[i];now=loss(z);out[i]+=now-previous;previous=now
    return out/6

out=[]
for case in ['V31_C60_warm_fine','V32_C60_warm_fine']:
 for method in ['T','N0','N10']:
    L=loadmat(SOURCE/case/(method+'.mat'),simplify_cells=True);P=prepare(L['C']);xs=L['Xs'];ts=L['Ts']
    values=[]
    for x in xs:
        V,ga,ch,_,_=actual_air(x);xf=x[19:25];e=np.array([V-xf[0],ga-xf[1],np.arctan2(np.sin(ch-xf[2]),np.cos(ch-xf[2]))])
        hr=calc(xf,np.zeros(3),P)[0][1];hc=calc(xf,e,P)[0][1];values.append(hr-hc)
    indices=sorted(set([len(xs)-1,int(np.argmax(values)),int(np.argmin(L['D'][:,15]))]))
    for i in indices:
        x=xs[i];xf=x[19:25];V,ga,ch,kap,rot=actual_air(x)
        e=np.array([V-xf[0],ga-xf[1],np.arctan2(np.sin(ch-xf[2]),np.cos(ch-xf[2]))])
        h0=calc(xf,np.zeros(3),P)[0];hc=calc(xf,e,P)[0];parts=attribution(xf,e,P)
        a=dict(case=case,method=method,index=i,time_s=float(ts[i]),is_last=i==len(xs)-1,is_max_loss=i==int(np.argmax(values)),
               reference_upper_N=h0[1],feedback_upper_N=hc[1],upper_margin_loss_N=h0[1]-hc[1],
               eV_mps=e[0],egamma_deg=e[1]*180/np.pi,echi_deg=e[2]*180/np.pi,
               shapley_loss_eV_N=parts[0],shapley_loss_egamma_N=parts[1],shapley_loss_echi_N=parts[2])
        assert abs(parts.sum()-(h0[1]-hc[1]))<1e-8
        out.append(a)
with (ROOT/'failure_attribution.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
(ROOT/'failure_attribution.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps([q for q in out if q['is_last'] and q['method']!='N10'],indent=2))
