"""Recheck W1-B interior counterexamples using exact outward intervals.
Run: python verify_w1b_points.py [output.json]
Only the declared nominal reference model is evaluated; no MATLAB/ODE run.
"""
import json
import sys
from pathlib import Path
from interval_certificate import I, SCALE, evaluate, ellipse, PV, PG, PC, q1, qg, qv, qc

CASES = {
    "A": ["0.56636307", "0.00239667", "0.23481519", "-0.03591428", "0"],
    "B": ["0.08392001", "-0.01059750", "-0.59489990", "-0.00555000", "0"],
}

def pack(x):
    return {"approximate": x.pair(), "lower_scaled_integer": x.lo,
            "upper_scaled_integer": x.hi, "scale": SCALE}

def main():
    result = {"scope": "nominal reference model only", "cases": {}}
    for name, coords in CASES.items():
        z = [I(x) for x in coords]
        bounds, root_valid = evaluate(z)
        vb = ellipse(z)
        if not root_valid or vb.hi >= I("0.20").lo:
            raise RuntimeError("Point is not verified inside the regular ellipsoid.")
        if name == "A":
            assert bounds["phiTmin"].hi < 0 and bounds["hTmin"].lo > 0
            assert all(x.lo > 0 for k, x in bounds.items() if k != "phiTmin")
        else:
            assert bounds["hTmin"].hi < 0 and bounds["phiTmin"].lo > 0
        result["cases"][name] = {
            "xi_decimal": coords, "VB": pack(vb),
            "Tff": pack(bounds["hTmin"] + 3),
            "hTmin": pack(bounds["hTmin"]), "phiTmin": pack(bounds["phiTmin"]),
            "all_checks": {k: pack(x) for k, x in bounds.items()},
        }
    # A rational lower bound for the generalized minimum eigenvalue.
    mu = I(".145"); a,b,c = PV; d,e,f = PG
    tests = {"speed_pivot": q1-mu*a,
             "speed_det": (q1-mu*a)*(qv-mu*c)-(mu*b).sq(),
             "climb_pivot": qg-mu*d,
             "climb_det": (qg-mu*d)*(qg-mu*f)-(mu*e).sq(),
             "course_pivot": qc-mu*PC}
    assert all(x.lo > 0 for x in tests.values())
    result["mu_lower_bound"] = "0.145"
    result["mu_verification"] = {k: pack(x) for k,x in tests.items()}
    out = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name("points_rechecked.json")
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("PASS: two interior counterexamples and mu_B > 0.145. Output:", out)

if __name__ == "__main__":
    main()
