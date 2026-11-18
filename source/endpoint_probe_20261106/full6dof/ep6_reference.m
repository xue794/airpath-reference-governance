function [r,v]=ep6_reference(t,C,O)
% Straight C2 deceleration, then constant-speed C2 course change.
% Position integrates exactly the same velocity law. During course change
% 24-point Gauss-Legendre quadrature is used; 48-point check is in the runner.
v0=C.R.Vb;vl=O.Vlow;t0=O.SpeedStart;dr=O.SpeedRise;
[s,I]=blend_integral(t,t0,dr);
speed=v0+(vl-v0)*s;
r=[v0*t+(vl-v0)*I;0;-C.R.h_ref];v=[speed;0;0];
if t<=O.TurnStart,return;end
turnEnd=O.TurnStart+O.TurnRise;
tpart=min(t,turnEnd);len=tpart-O.TurnStart;
nodes=O.Nodes;weights=O.Weights;
z=O.TurnStart+.5*len*(nodes+1);
q=(z-O.TurnStart)/O.TurnRise;
psi=O.Course*(10*q.^3-15*q.^4+6*q.^5);
delta=.5*len*vl*[sum(weights.*(cos(psi)-1));sum(weights.*sin(psi));0];
r=r+delta;
if t>turnEnd
    r=r+(t-turnEnd)*vl*[cos(O.Course)-1;sin(O.Course);0];
end
q=min(1,(t-O.TurnStart)/O.TurnRise);
psi=O.Course*(10*q^3-15*q^4+6*q^5);
v=vl*[cos(psi);sin(psi);0];
end
function [s,I]=blend_integral(t,t0,d)
if t<=t0,s=0;I=0;elseif t>=t0+d,s=1;I=t-t0-d/2;else
q=(t-t0)/d;s=10*q^3-15*q^4+6*q^5;I=d*(2.5*q^4-3*q^5+q^6);end
end
