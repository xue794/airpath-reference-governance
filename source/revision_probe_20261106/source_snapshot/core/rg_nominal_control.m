function u0 = rg_nominal_control(x,cc,R)
%RG_NOMINAL_CONTROL Unconstrained second-order reference filter.
eta=[x(1)-cc(1);
     x(2)-cc(2);
     rg_angle_error(x(3),cc(3))];
nu=x(4:6);
u0=-2*R.Zeta*R.Omega*nu-(R.Omega^2)*eta;
end
