function Cc=viB_command_attitude(gamma,chi,mu,alpha)
%VIB_COMMAND_ATTITUDE Command body-to-inertial DCM.
%
% C_IB,c = Rz(chi) Ry(gamma) Rx(mu) Ry(alpha)
Cc=viB_Rz(chi)*viB_Ry(gamma)*viB_Rx(mu)*viB_Ry(alpha);
end
