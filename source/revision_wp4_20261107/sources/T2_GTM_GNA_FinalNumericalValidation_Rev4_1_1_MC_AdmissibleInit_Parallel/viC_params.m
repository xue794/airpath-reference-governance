function F = viC_params()
%VIC_PARAMS Section VI-C SO(3) command-filter module parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% FINAL SYNCHRONIZED parameter values for the Rev.1 regression freeze.
% The primary certified requirement is that the realized command-filter
% attitude error remains inside theta_f <= 0.30 deg.

B=viB_params();

F.revision='Section VI-C SO(3) Command Filter -- Final Synchronized Rev.1';

F.kRf=256.0;
F.KomegaF=7.0*eye(3);
F.cF=0.50;

F.thetaFMax=B.thetaFMax;
F.PsiFMax=1-cos(F.thetaFMax);

% Quaternion norm stabilization vanishes identically on S^3.
F.quatNormGain=1.0;

% Constant-command verification.
F.constT=5.0;
F.constOffset=deg2rad(0.25);
F.constOmega0=[0.020;-0.015;0.010];

% Persistent pure-yaw command used to demonstrate nonzero moving-reference lag.
F.spinRate=0.030;          % rad/s
F.spinT=80.0;
F.spinSteadyT=60.0;

% Numerical checks.
F.fdDt=1e-7;
F.Nkin=300;
F.Nlyap=5000;
F.RelTol=1e-9;
F.AbsTol=1e-11;
F.filterMaxStep=0.02;
end
