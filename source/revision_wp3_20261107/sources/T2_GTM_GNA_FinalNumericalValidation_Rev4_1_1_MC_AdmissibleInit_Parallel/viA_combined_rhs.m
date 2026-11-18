function zdot = viA_combined_rhs(t,z,mode,R,A,E)
%VIA_COMBINED_RHS Section V governor + Section VI-A ideal air-path subsystem.
%
% State:
%   z = [xf(6); y(3)]
%
% No attitude dynamics, sideslip dynamics, actuator dynamics, or real T-2
% force/moment realization are included in VI-A.

xf=z(1:6);
y=z(7:9);

S=viA_case_signal(t,mode,E);
[vc,~]=rg_raw_test_command(t,S.refMode,R);
[cc,~]=rg_raw_command_chart(vc,R);

u0=rg_nominal_control(xf,cc,R);
O=rg_scalar_governor(xf,u0,R);
xfDot=[xf(4:6);O.u];

yDot=viA_ideal_path_rhs(y,xf,S.wdot,R,A);

zdot=[xfDot;yDot];
end
