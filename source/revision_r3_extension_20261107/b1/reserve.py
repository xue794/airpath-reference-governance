"""B1 one-box analytic reserve certificate for the frozen binary64 nominal model.

Rational operations are exact. Square roots use checked outward rational bounds.
The B2 subdivision, grid search and coverage audit are separate.
"""
from fractions import Fraction as Q
from math import isqrt

def binary(x):
    return Q.from_float(float(x))

def sqrt_upper(x):
    unit=10**45
    n=isqrt(x.numerator*unit**2//x.denominator)
    upper=Q(n,unit)
    if upper**2<x:
        upper+=Q(1,unit)
    assert upper**2>=x and (upper==0 or (upper-Q(1,unit))**2<x)
    return upper

def decimal_bound(q,places,upper):
    scale=10**places
    i=-(-q*scale//1) if upper else q*scale//1
    sign="-" if i<0 else ""
    a,b=divmod(abs(i),scale)
    return f"{sign}{a}.{b:0{places}d}"

def inverse_diagonal(P):
    n=len(P)
    assert all(P[i][j]==P[j][i] for i in range(n) for j in range(n))
    aug=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(P)]
    for i in range(n):
        # For a symmetric matrix the successive positive Schur pivots certify SPD.
        pivot=aug[i][i]
        assert pivot>0
        aug[i]=[v/pivot for v in aug[i]]
        for j in range(n):
            if i!=j:
                v=aug[j][i]
                aug[j]=[a-v*b for a,b in zip(aug[j],aug[i])]
    assert all(aug[i][j]==(i==j) for i in range(n) for j in range(n))
    return [aug[i][n+i] for i in range(n)]

def coeff_at(coeff,a):
    return sum(v*a**i for i,v in enumerate(coeff))

def polynomial_range(coeff,ahi):
    """Exact power-by-power enclosure on [0,ahi]."""
    low=high=coeff[0]
    for i,c in enumerate(coeff[1:],1):
        end=c*ahi**i
        low+=min(0,end); high+=max(0,end)
    return low,high

def uniform(C):
    R=C["R"]; A=C["A"]; P=R["P"]
    mass=binary(P["m"]); grav=binary(R["g"]); K=binary(R["rho"])*binary(P["S"])/2
    v0=binary(R["Vb"]); alo=binary(R["alphaLo"]); ahi=binary(R["alphaHi"])
    reserve=binary(R["DeltaN"]); cap=Q(1)
    assert alo==0 and 0<ahi<Q(1,5) and R["gammaB"]==0 and reserve==10 and v0==binary(38.1)
    # Unify the stipulated decimal error scale and its binary64 implementation.
    scale=min(Q(381,10),v0)
    rho=binary(R["rhoB"]); assert rho==binary(.1)
    pb=[[binary(x) for x in row] for row in R["PB"]]
    h=[sqrt_upper(rho*d) for d in inverse_diagonal(pb)]
    hv,hg,hnv,hng,hnc=h
    vflo=v0-hv; vfhi=v0+hv; vlo=vflo-cap; vhi=vfhi+cap
    max_g=hg+cap/scale
    assert vlo>0 and max_g<1  # |gamma|<1 rad implies cos(gamma)>0.
    cl=[binary(P["theta"][n]) for n in [15,16,20,21,22]]
    cd=[binary(P["theta"][n]) for n in [0,1,4,7,9]]
    kv,kg,kc=[binary(A[k]) for k in ("kV","kGamma","kChi")]
    assert min(kv,kg,kc)>0
    # Use cos(a)>=1-a^2/2, 0<=a<=tan(a)<=a/(1-a^2/2).
    cos_min=1-ahi**2/2; tan_hi=ahi/cos_min
    CLhi=coeff_at(cl,ahi); CDhi=coeff_at(cd,ahi)
    assert CLhi>0 and CDhi>0 and cl[0]>0
    upper_B=CLhi+CDhi*tan_hi
    # Increment coefficients per unit joint error radius.
    dX=mass*(kv+grav/scale)
    dAg=hng+(vhi*kg+grav*max_g)/scale
    dVcos=1+vfhi*max_g/scale
    dAc=hnc*dVcos+vhi*kc/scale
    dN=mass*sqrt_upper(dAg**2+dAc**2)
    dQ=K*(2*vfhi+cap)
    Llo=dN+dQ*abs(cl[0]); Lhi=dN+dQ*upper_B+dX*tan_hi
    # min/max raw reference endpoint gaps, independently of any old certificate.
    Xfabs=mass*(hnv+grav*hg)
    agmin=grav*(1-hg**2/2)-vfhi*hng; assert agmin>0
    Nflo=mass*agmin
    Nfhi=mass*sqrt_upper((grav+vfhi*hng)**2+(vfhi*hnc)**2)
    gflo=Nflo-K*vfhi**2*cl[0]
    gfhi=K*vflo**2*(CLhi+CDhi*ahi)-Xfabs*tan_hi-Nfhi
    assert min(gflo,gfhi)>reserve
    # dM/da for every candidate angle over the full 1 m/s auxiliary domain.
    dcl_low,_=polynomial_range([i*c for i,c in enumerate(cl) if i],ahi)
    cdr=polynomial_range(cd,ahi); dcdr=polynomial_range([i*c for i,c in enumerate(cd) if i],ahi)
    maxcd=max(map(abs,cdr)); maxdcd=max(map(abs,dcdr))
    Xabs=mass*(hnv+kv*cap+grav*max_g)
    assert dcl_low>0
    Mprime=K*vlo**2*dcl_low-K*vhi**2*maxdcd*tan_hi-(K*vhi**2*maxcd+Xabs)/cos_min**2
    assert Mprime>0
    threshold=min(cap,reserve/max(Llo,Lhi))
    radius=Q(-(-1000*threshold//1)-1,1000)  # largest 0.001-grid point strictly below threshold
    assert radius>0 and max(Llo,Lhi)*radius<reserve
    return dict(projection=h,Llower=Llo,Lupper=Lhi,D_X=dX,D_Agamma=dAg,D_Achi=dAc,D_N=dN,
        D_q=dQ,B_upper=upper_B,V_lower=vlo,gamma_abs=max_g,Mprime_lower=Mprime,
        raw_ref_lower=gflo,raw_ref_upper=gfhi,threshold=threshold,radius=radius,
        scale=scale,reserve=reserve)

def display(cert):
    upper=("Llower","Lupper","D_X","D_Agamma","D_Achi","D_N","D_q","B_upper","gamma_abs")
    d={k:decimal_bound(v,12,k in upper) for k,v in cert.items() if k!="projection"}
    d["projection"]=[decimal_bound(v,18,True) for v in cert["projection"]]
    d["examples"]=[]
    for eps in [Q(0),Q(1,8),cert["radius"],Q(1)]:
        lo,hi=[eps*cert[k] for k in ["Llower","Lupper"]]
        d["examples"].append(dict(epsilon_mps=decimal_bound(eps,6,False),
            angle_half_rad_upper=decimal_bound(eps/cert["scale"],12,True),
            U_lower_N=decimal_bound(lo,9,True),U_upper_N=decimal_bound(hi,9,True),
            ten_remaining_lower_N=decimal_bound(10-lo,9,False),ten_remaining_upper_N=decimal_bound(10-hi,9,False),
            KB_remaining_lower_N=decimal_bound(cert["raw_ref_lower"]-lo,9,False),
            KB_remaining_upper_N=decimal_bound(cert["raw_ref_upper"]-hi,9,False)))
    return d
