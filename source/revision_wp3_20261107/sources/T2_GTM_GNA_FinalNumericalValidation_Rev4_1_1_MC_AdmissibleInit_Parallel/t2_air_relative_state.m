function A = t2_air_relative_state(t,x,U,isQuat)
%T2_AIR_RELATIVE_STATE Air-relative diagnostic quantities for Euler/quaternion states.
if nargin<4, isQuat=true; end
if isQuat
    CIB=quat_to_dcm(x(7:10)); vgB=x(4:6); omega=x(11:13);
else
    CIB=euler321_dcm(x(7),x(8),x(9)); vgB=x(4:6); omega=x(10:12);
end
[wI,wDotI]=t2_wind_eval(t,U);
vaB=vgB-CIB.'*wI; V=norm(vaB);
alpha=atan2(vaB(3),vaB(1)); beta=asin(max(-1,min(1,vaB(2)/V)));
A=struct('CIB',CIB,'vgB',vgB,'vaB',vaB,'V',V,'alpha',alpha,'beta',beta, ...
         'omega',omega,'wI',wI,'wDotI',wDotI);
end
