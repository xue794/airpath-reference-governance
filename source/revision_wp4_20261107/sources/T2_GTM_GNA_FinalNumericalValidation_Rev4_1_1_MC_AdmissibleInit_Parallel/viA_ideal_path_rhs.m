function [ydot,C] = viA_ideal_path_rhs(y,xf,wdot,R,A)
%VIA_IDEAL_PATH_RHS Exact ideal coordinated air-path dynamics.
%
% Implements
%   Vdot = -g sin(gamma) + (Fparallel-D0)/m - t'*wdot
%   gammadot = -g cos(gamma)/V + N cos(mu)/(mV) - nGamma'*wdot/V
%   chidot = N sin(mu)/(mV cos(gamma)) - nChi'*wdot/(V cos(gamma))

C=viA_feedback_command(y,xf,R,A);

V=y(1); gamma=y(2); chi=y(3);
B=viA_airpath_basis(gamma,chi);

Vdot=-R.g*sin(gamma)+(C.Fparallel-C.D)/R.P.m-B.t.'*wdot;

gammaDot=-R.g*cos(gamma)/V + ...
         C.N*cos(C.mu)/(R.P.m*V) - ...
         B.nGamma.'*wdot/V;

chiDot=C.N*sin(C.mu)/(R.P.m*V*cos(gamma)) - ...
       B.nChi.'*wdot/(V*cos(gamma));

ydot=[Vdot;gammaDot;chiDot];
end
