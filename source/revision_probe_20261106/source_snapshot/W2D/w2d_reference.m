function [rd,rdDot,M]=w2d_reference(t,mode,C,O)
%W2D_REFERENCE Smooth ground-position reference whose zero-error still-wind
% velocity reaches the predeclared W2-B raw air-path witness.
if nargin<4||isempty(O),O=w2d_options();end
K=w2d_case_info(mode,C,O);h0=C.R.h_ref;
Vb=C.R.Vb;v0=[Vb;0;0];
Vt=K.TargetAirPath(1);gt=K.TargetAirPath(2);ct=K.TargetAirPath(3);
vt=Vt*[cos(gt)*cos(ct);cos(gt)*sin(ct);-sin(gt)];
[s,I]=w2d_smooth_pulse(t,K.RampStart_s,K.RampDuration_s,K.PulseHold_s);
rd0=[0;0;-h0];rd=rd0+v0*t+(vt-v0)*I;rdDot=v0+(vt-v0)*s;
M=struct('Blend',s,'BlendIntegral_s',I,'InitialVelocity',v0,'TargetVelocity',vt,'TargetAirPath',K.TargetAirPath);
end
