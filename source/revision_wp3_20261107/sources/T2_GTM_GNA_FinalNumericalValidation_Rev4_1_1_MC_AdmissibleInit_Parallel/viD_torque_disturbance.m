function d=viD_torque_disturbance(t,mode,Dp)
%VID_TORQUE_DISTURBANCE Deterministic bounded rotational residual [N m].

if strcmpi(mode,'disturbed')
    d=Dp.distAmp.*sin(Dp.distOmega*t+Dp.distPhase);
else
    d=zeros(3,1);
end
end
