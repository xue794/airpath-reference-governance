function zDot = wind_estimator_augmented_rhs(t,z,U,P,E,addNoise)
%WIND_ESTIMATOR_AUGMENTED_RHS Quaternion T-2 plant + wind estimator.
% z = [xPlant(13); what(3)].
x=z(1:13); what=z(14:16);
xd=t2_gna_6dof_quat(t,x,U,P);
CIB=quat_to_dcm(x(7:10));
[wI,~]=t2_wind_eval(t,U);
M=wind_estimator_measurement_from_truth(CIB,x(4:6),wI,t,E,addNoise);
whatDot=E.Lw*(M.yw-what);
zDot=[xd;whatDot];
end
