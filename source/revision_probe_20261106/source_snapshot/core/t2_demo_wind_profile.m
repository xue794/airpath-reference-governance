function z = t2_demo_wind_profile(t)
%T2_DEMO_WIND_PROFILE Smooth time-varying inertial wind for Rev.1 demos.
% Output z=[wI;wDotI], both in NED inertial axes.
w=[4+2*sin(0.4*t); -3+1.5*cos(0.7*t); 0.5*sin(0.3*t)];
wd=[0.8*cos(0.4*t); -1.05*sin(0.7*t); 0.15*cos(0.3*t)];
z=[w;wd];
end
