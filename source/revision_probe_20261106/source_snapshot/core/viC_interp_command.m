function Cc=viC_interp_command(t,Cmd)
%VIC_INTERP_COMMAND Interpolate the already-verified VI-A command trajectory.
%
% The command filter does not feed back into VI-A, so the upstream trajectory
% can be generated once and interpolated during the filter integration.

gamma=interp1(Cmd.t,Cmd.gamma,t,'pchip');
chi=interp1(Cmd.t,Cmd.chi,t,'pchip');
mu=interp1(Cmd.t,Cmd.mu,t,'pchip');
alpha=interp1(Cmd.t,Cmd.alpha,t,'pchip');

Cc=viB_command_attitude(gamma,chi,mu,alpha);
end
