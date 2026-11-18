"""Independent arithmetic smoke tests; not a formal software proof.
Standard-library rational tests always run. Optional mpmath tests use 80 digits.
"""
from fractions import Fraction
from random import Random
import json
from pathlib import Path
from interval_certificate import I, SCALE, sin, cos, PI

def contains(box, x):
    return Fraction(box.lo, SCALE) <= x <= Fraction(box.hi, SCALE)

def main():
    rng=Random(32026); count=0
    for _ in range(100):
        aa=sorted([Fraction(rng.randint(-900,900),1000) for _ in range(2)])
        bb=sorted([Fraction(rng.randint(100,900),1000) for _ in range(2)])
        a=I(str(aa[0]),str(aa[1])); b=I(str(bb[0]),str(bb[1]))
        for box, fn in [(a+b,lambda x,y:x+y),(a-b,lambda x,y:x-y),
                        (a*b,lambda x,y:x*y),(a/b,lambda x,y:x/y)]:
            assert all(contains(box,fn(x,y)) for x in aa for y in bb)
            count+=1
    result={"rational_corner_tests":count,"status":"PASS"}
    try:
        import mpmath as mp
    except ImportError:
        result["transcendental_independent_tests"]="SKIPPED: mpmath not installed"
    else:
        mp.mp.dps=80
        inside=lambda v,x: mp.mpf(v.lo)/SCALE <= x <= mp.mpf(v.hi)/SCALE
        assert inside(PI,mp.pi)
        for _ in range(100):
            x=str(Fraction(rng.randint(-900,900),1000)); q=I(x); exact=mp.mpf(q.lo)/SCALE
            assert inside(sin(q),mp.sin(exact))
            assert inside(cos(q),mp.cos(exact))
        result.update({"transcendental_independent_tests":200,"pi_enclosure":True})
    Path(__file__).with_name("selftest_result.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
