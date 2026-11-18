function O=vii_actuator_derivative(delta,T,deltaCmd,Tcmd,P7)
%VII_ACTUATOR_DERIVATIVE Bounded first-order surface/propulsor interface.

deltaCmd=min(max(deltaCmd,P7.deltaMin),P7.deltaMax);
Tcmd=min(max(Tcmd,P7.Tmin),P7.Tmax);

dDelta=(deltaCmd-delta)./P7.tauDelta;
dDelta=min(max(dDelta,-P7.deltaRateMax),P7.deltaRateMax);

% Forward-invariance projection at the hard position box.
for i=1:3
    if delta(i)<=P7.deltaMin(i) && dDelta(i)<0
        dDelta(i)=0;
    elseif delta(i)>=P7.deltaMax(i) && dDelta(i)>0
        dDelta(i)=0;
    end
end

dT=(Tcmd-T)/P7.tauT;
dT=min(max(dT,-P7.TRateMax),P7.TRateMax);

if T<=P7.Tmin && dT<0, dT=0; end
if T>=P7.Tmax && dT>0, dT=0; end

O.deltaDot=dDelta;
O.TDot=dT;
O.deltaCmd=deltaCmd;
O.Tcmd=Tcmd;
end
