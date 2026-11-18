"""Frozen breadth-first grid, eight refinement steps and final rerun."""
import csv,gzip,hashlib,io,json
from collections import Counter,deque
from fractions import Fraction as F
from pathlib import Path
from boxes import Kernel,S,H,ceil,MODEL,B1

GRID=list(map(F,["0","1/64","1/16","1/8","1/4","1/2","3/4","1"]))
BUDGET=100000
OUT=H/"trees"
def tag(e):return f"{e.numerator}_{e.denominator}"
def writej(p,data):p.write_text(json.dumps(data,indent=2,allow_nan=False),encoding="utf8")
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def search(e,sub=False,repeat=False):
    stem=("sub" if sub else "full")+"_"+tag(e)+("_repeat" if repeat else "")
    k=Kernel(e,sub);todo=deque([1]);counts=Counter();generated=1
    path=OUT/(stem+".csv.gz");mx=[0,0];worst=[None,None];visited=0
    with path.open("wb") as fd,gzip.GzipFile(filename="",fileobj=fd,mode="wb",mtime=0,compresslevel=3) as gz,io.TextIOWrapper(gz,encoding="ascii",newline="") as f:
        writer=csv.writer(f,lineterminator="\n")
        writer.writerow(["heap_id","state","upper_lower_pN","upper_upper_pN","domain_lower_1e15"])
        while todo:
            node=todo.popleft();kind,l,u,dl=k.node(node);visited+=1
            if kind=="U" and generated+2<=BUDGET:
                kind="S";todo.extend([2*node,2*node+1]);generated+=2
            if kind in ("A","U"):
                for d,v in enumerate([l,u]):
                    if worst[d] is None or v>mx[d]:mx[d]=v;worst[d]=node
            counts[kind]+=1
            writer.writerow([node,kind,ceil(l*10**12,S),ceil(u*10**12,S),dl*10**15//S])
    assert visited==generated==2*counts["S"]+1 and counts["P"]+counts["A"]+counts["U"]==counts["S"]+1
    row=dict(domain="subset" if sub else "full",epsilon=str(e),eps=float(e),tree=path.name,sha256=digest(path),
        visited=visited,pruned=counts["P"],accepted=counts["A"],split=counts["S"],unresolved=counts["U"],
        max_lower_pN=ceil(mx[0]*10**12,S),max_upper_pN=ceil(mx[1]*10**12,S),worst_ids=worst,
        passed=counts["U"]==0,repeat=repeat)
    print(stem,visited,"unresolved",counts["U"],"U",row["max_lower_pN"]/1e12,row["max_upper_pN"]/1e12,flush=True)
    k.axis.cache_clear()
    return row

def run():
    OUT.mkdir(exist_ok=True); rows=[]
    for sub in [False,True]:
        grid=[search(e,sub) for e in GRID];rows.extend(grid)
        passed=[F(r["epsilon"]) for r in grid if r["passed"]]
        assert F(0) in passed
        lo=max(passed)
        if lo==0 and not sub:continue
        if 0<lo<1:
            hi=GRID[GRID.index(lo)+1]
            for _ in range(8):
                mid=(hi+lo)/2;r=search(mid,sub);rows.append(r)
                if r["passed"]:lo=mid
                else:hi=mid
        last=search(lo,sub,True);rows.append(last)
        old=next(row for row in rows if row["domain"]==last["domain"] and F(row["epsilon"])==lo)
        assert last["sha256"]==old["sha256"] and last["passed"]
        break
    writej(H/"search.json",dict(rows=rows,certified_halfwidth=str(lo),domain=last["domain"],b1_hash=digest(B1),
        regularity_source="B1 full box, all eps <= 1; binary64 nominal kernel",
        inherited_regular={k:MODEL["analytic"][k] for k in ("Mprime_lower","V_lower","gamma_abs","raw_ref_lower","raw_ref_upper")}))
if __name__=="__main__":run()
