function [fX,X0]=viii_reduced_vector_field(t,X,mode,C,P8)
%VIII_REDUCED_VECTOR_FIELD Directional numerical derivative of X_red.
%
% The complete 36-state RHS is the source of truth.  A central directional
% derivative of the reduced-coordinate map is used so command derivatives,
% actuator saturation, wind kinematics, and rigid-body couplings are all
% retained without introducing a separate surrogate model.

if nargin<5 || isempty(P8), P8=viii_params(); end

[x,~]=viii_lift_reduced(t,X,mode,C,P8);
[X0,~]=viii_reduced_from_full(t,x,mode,C,P8);
rhs=cl6_rhs(t,x,mode,C);

h=P8.flowDt;

xp=x+h*rhs;
xm=x-h*rhs;

% Normalize the two quaternions after the Euler directional step.
xp(7:10)=quat_normalize(xp(7:10));
xm(7:10)=quat_normalize(xm(7:10));
xp(26:29)=quat_normalize(xp(26:29));
xm(26:29)=quat_normalize(xm(26:29));

[Xp,~]=viii_reduced_from_full(t+h,xp,mode,C,P8);
[Xm,~]=viii_reduced_from_full(t-h,xm,mode,C,P8);

fX=(Xp-Xm)/(2*h);
end
