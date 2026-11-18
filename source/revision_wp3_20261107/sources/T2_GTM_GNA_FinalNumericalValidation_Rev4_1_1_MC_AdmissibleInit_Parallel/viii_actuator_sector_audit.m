function A=viii_actuator_sector_audit(C,P8)
%VIII_ACTUATOR_SECTOR_AUDIT Exact sector test for the implemented rate-limited actuators.
%
% With command held fixed, e=act-cmd satisfies
%   e_i*phi_i(e_i) >= kappa_i e_i^2
% on the full physical box, where phi is the sign-preserving rate-limited
% first-order restoring map.  The hard-box projection can only remove an
% outward-pointing velocity and therefore does not weaken this property
% for commands inside the same box.

N=P8.Nsector;
minRatio=inf;
maxViolation=0;

for i=1:3
    emax=C.P7.deltaMax(i)-C.P7.deltaMin(i);
    eg=linspace(-emax,emax,N);
    for k=1:N
        e=eg(k);
        if abs(e)<1e-14, continue; end
        phi=min(max(e/C.P7.tauDelta(i),-C.P7.deltaRateMax(i)), ...
                C.P7.deltaRateMax(i));
        ratio=e*phi/(e^2);
        minRatio=min(minRatio,ratio);
        maxViolation=max(maxViolation,P8.kappaDelta(i)-ratio);
    end
end

emaxT=C.P7.Tmax-C.P7.Tmin;
eg=linspace(-emaxT,emaxT,N);
for k=1:N
    e=eg(k);
    if abs(e)<1e-14, continue; end
    phi=min(max(e/C.P7.tauT,-C.P7.TRateMax),C.P7.TRateMax);
    ratio=e*phi/(e^2);
    minRatio=min(minRatio,ratio);
    maxViolation=max(maxViolation,P8.kappaT-ratio);
end

A.minRatio=minRatio;
A.maxViolation=maxViolation;
A.kappaDelta=P8.kappaDelta;
A.kappaT=P8.kappaT;
A.kappaA=P8.kappaA;
end
