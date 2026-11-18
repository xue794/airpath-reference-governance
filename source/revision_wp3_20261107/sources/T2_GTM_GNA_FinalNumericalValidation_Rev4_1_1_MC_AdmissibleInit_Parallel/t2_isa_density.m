function rho = t2_isa_density(h)
%T2_ISA_DENSITY  ISA tropospheric density for h in meters.
T0=288.15; p0=101325; L=0.0065; R=287.05287; g=9.80665;
if h < 0
    h = 0;
end
if h > 11000
    error('t2_isa_density:h','This simple routine is limited to h <= 11 km.');
end
T = T0-L*h;
p = p0*(T/T0)^(g/(R*L));
rho = p/(R*T);
end
