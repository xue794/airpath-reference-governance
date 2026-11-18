"""Independently check topology and every leaf of each claimed passing tree."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from collections import Counter
import csv,gzip,json
from independent import Check,iv,q,fraction
from search import H,OUT,digest,writej

def batch(arg):
    e,sub,rows=arg;c=Check(e,sub)
    for n,kind,ul,uu in rows:
        z=c.bounds(n)
        if kind=="P":
            assert c.domain(z)>c.r,(e,n)
        else:
            losses=c.coeff_losses(z)
            for i,(v,bound) in enumerate(zip(losses,(ul,uu))):
                u=min(c.caps[i].b,max(q(0),v.b))
                assert u<10 and u<=q(F(bound,10**12)).b,(e,n,i,str(u),bound)
    return len(rows)

def review():
    search=json.loads((H/"search.json").read_bytes());out=[];jobs=[]
    for row in search["rows"]:
        p=OUT/row["tree"];assert digest(p)==row["sha256"]
        c=Check(row["epsilon"],row["domain"]=="subset");needed={1};cnt=Counter();ums=[0,0];worst=[]
        with gzip.open(p,"rt") as f:
            chunks=[]
            for r in csv.DictReader(f):
                n=int(r["heap_id"]);s=r["state"];needed.remove(n);cnt[s]+=1
                lu,uu=int(r["upper_lower_pN"]),int(r["upper_upper_pN"])
                if s=="S":needed.update([2*n,2*n+1])
                if s in ("A","U"):ums=[max(ums[0],lu),max(ums[1],uu)]
                if row["passed"] and s in ("A","P") and not row["repeat"]:
                    chunks.append((n,s,lu,uu))
                    if len(chunks)==512:
                        jobs.append((c.eps,row["domain"]=="subset",chunks));chunks=[]
                if n in row["worst_ids"]:
                    bounds=c.box(n);point=[q((l+h)/2) for l,h in bounds]
                    vals=c.direct(point);enc=c.coeff_losses(c.bounds(n))
                    assert all(a.a>=b.a-q("1e-48") and a.b<=b.b+q("1e-48") for a,b in zip(vals,enc))
                    worst.append(dict(heap_id=n,center_loss_N=[float(a.mid) for a in vals],
                        center_membership_ratio=float((c.domain(point)/c.r).mid),
                        ledger_pN=[lu,uu],center_coordinates=[str((l+h)/2) for l,h in bounds]))
            if chunks:jobs.append((c.eps,row["domain"]=="subset",chunks))
        assert not needed and sum(cnt.values())==row["visited"]
        assert ums==[row["max_lower_pN"],row["max_upper_pN"]]
        if row["repeat"]:
            assert any(x["sha256"]==row["sha256"] and not x["repeat"] for x in search["rows"])
        out.append(dict(tree=p.name,topology=True,counts=dict(cnt),worst_center_regressions=worst))
    with ProcessPoolExecutor(12) as pool: checked=sum(pool.map(batch,jobs,chunksize=1))
    assert checked==sum(r["accepted"]+r["pruned"] for r in search["rows"] if r["passed"] and not r["repeat"])
    writej(H/"tree_review.json",dict(checked_distinct_passing_tree_leaves=checked,repeat_verified_byte_identical=True,trees=out))
    print("Independent leaf and topology checks passed",checked)
if __name__=="__main__":review()
