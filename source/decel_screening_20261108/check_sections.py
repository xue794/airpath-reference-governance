"""Intersect the already frozen cross-section points with the feedback window."""
from scipy.io import loadmat
from freeze import lock,OLD,RESULT
from screen_feedback import ForceCheck,guide,pt,BOXES,csvout

def run():
    lock()
    model=loadmat(OLD/'N10.mat',simplify_cells=True)
    probe=loadmat(RESULT/'affine_export.mat',simplify_cells=True)
    engine=ForceCheck(model['C']); rows=[]
    for name,s in zip(probe['labels'],probe['sets']):
        if not name.startswith('section'):continue
        d=int(name[7:9])
        for i,t in enumerate(s['t']):
            if not 14<=t<=20:continue
            center=guide(float(t),d);offs=s['offset'][i]
            x=tuple(c+pt(e) for c,e in zip(center,[*offs[:2],0,*offs[2:]]))
            for b in BOXES:rows.append(dict(case=name,sample_index=i,time_s=t,**engine.examine(x,b)))
    csvout('feedback_section_probe.csv',rows)

if __name__=='__main__':run()
