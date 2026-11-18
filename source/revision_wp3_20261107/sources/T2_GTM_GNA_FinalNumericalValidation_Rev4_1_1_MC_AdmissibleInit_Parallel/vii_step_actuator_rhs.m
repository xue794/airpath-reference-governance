function xdot=vii_step_actuator_rhs(~,x,deltaCmd,Tcmd,P7)
%VII_STEP_ACTUATOR_RHS Standalone actuator box/rate stress test.
O=vii_actuator_derivative(x(1:3),x(4),deltaCmd,Tcmd,P7);
xdot=[O.deltaDot;O.TDot];
end
