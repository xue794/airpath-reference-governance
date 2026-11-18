function S=viB_synthetic_relative_rotations(t,B)
%VIB_SYNTHETIC_RELATIVE_ROTATIONS Smooth bounded Qf/QR test trajectories.

af=[1.0;
    0.40*sin(0.11*t);
    0.30*cos(0.07*t)];
af=af/norm(af);

ar=[0.30*cos(0.09*t);
    0.80;
    0.40*sin(0.05*t+0.20)];
ar=ar/norm(ar);

angF=B.syntheticScaleF*B.thetaFMax*sin(0.17*t);
angR=B.syntheticScaleR*B.thetaRMax*cos(0.13*t+0.20);

S.Qf=rotvec_to_dcm(af*angF);
S.QR=rotvec_to_dcm(ar*angR);
S.thetaF=abs(angF);
S.thetaR=abs(angR);
end
