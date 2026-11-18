function M = wind_estimator_measurement_from_truth(CIB,vgB,wI,t,E,addNoise)
%WIND_ESTIMATOR_MEASUREMENT_FROM_TRUTH Build pseudo-wind measurement from plant truth.
% CIB maps B->I. vgB is ground-relative velocity in B. wI is inertial wind.
if nargin<6, addNoise=false; end
vgI=CIB*vgB;
vaB=vgB-CIB.'*wI;
z=vaB_to_airdata(vaB);
if addNoise
    N=wind_estimator_test_noise(t,E);
    vgIm=vgI+N.ng;
    % Measurement attitude error: C_m = C * exp(hat(dtheta)).
    CIBm=CIB*rotvec_to_dcm(N.dtheta);
    zm=z+N.nz;
else
    N=struct('ng',zeros(3,1),'dtheta',zeros(3,1),'nz',zeros(3,1));
    vgIm=vgI; CIBm=CIB; zm=z;
end
[yw,nw]=wind_pseudomeasurement(vgIm,CIBm,zm,wI);
M=struct('vgI',vgI,'vaB',vaB,'z',z,'vgIm',vgIm,'CIBm',CIBm,'zm',zm, ...
         'yw',yw,'nw',nw,'noise',N);
end
