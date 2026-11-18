function e = rg_angle_error(a,b)
%RG_ANGLE_ERROR Local S1 error a-b in (-pi,pi].
e = atan2(sin(a-b),cos(a-b));
end
