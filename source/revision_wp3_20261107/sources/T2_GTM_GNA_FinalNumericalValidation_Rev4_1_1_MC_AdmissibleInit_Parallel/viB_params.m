function B = viB_params()
%VIB_PARAMS Section VI-B command-attitude geometry / AoA-sideslip transfer.
% MATLAB R2018a compatible; base MATLAB only.
%
% This is a LOCAL certified attitude tube.  Its size is selected so that the
% Section V alpha_c interval [0,9] deg plus the geometric alpha deviation
% remains strictly inside the physical [-2,10] deg envelope.

R=rg_params();

B.revision='Section VI-B Command Attitude Geometry + AoA/Sideslip Transfer Rev.0';

B.thetaFMax=deg2rad(0.30);
B.thetaRMax=deg2rad(0.40);
B.thetaTube=B.thetaFMax+B.thetaRMax;
B.thetaBar=deg2rad(5.0);

B.betaBar=B.thetaTube;
B.deltaAlphaMax=B.thetaTube/cos(B.thetaBar);

B.alphaCBar=R.alphaHi;
B.kappaMin=cos(B.alphaCBar+B.deltaAlphaMax)*cos(B.betaBar);

% Numerical audit sizes.
B.Nmc=10000;
B.Nder=500;
B.fdDt=1e-7;

% Synthetic dynamic relative-attitude errors use the full local tube amplitudes.
B.syntheticScaleF=1.0;
B.syntheticScaleR=1.0;
end
