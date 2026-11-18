function xdot = rg_governor_rhs(t,x,mode,R)
%RG_GOVERNOR_RHS Second-order certified reference governor.
[vc,~]=rg_raw_test_command(t,mode,R);
[cc,~]=rg_raw_command_chart(vc,R);
u0=rg_nominal_control(x,cc,R);
O=rg_scalar_governor(x,u0,R);
xdot=[x(4:6);O.u];
end
