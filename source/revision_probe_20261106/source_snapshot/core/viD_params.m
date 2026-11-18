function Dp=viD_params()
%VID_PARAMS Section VI-D ideal rotational PHS tracking parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% This revision realizes the virtual rotational PHS torque port directly.
% Surface allocation, actuator dynamics, tau_0, S_omega,a^nom compensation,
% and propulsion moments remain Section VII tasks.

B=viB_params();

Dp.revision='Section VI-D PHS-Based Rotational Attitude/Rate Controller Rev.0';

% Rotational tracking gains.
Dp.kR=80.0;                 % N m / rad-equivalent geometric error
Dp.KomegaC=6.0*eye(3);      % N m s / rad
Dp.cR=0.25;                 % local storage cross-term coefficient

% The dynamic VI-B/VI-C certified rotational tube.
Dp.thetaRMax=B.thetaRMax;   % 0.40 deg

% Initial errors used by module tests.
Dp.constTheta0=deg2rad(0.30);
Dp.constEOmega0=[0.020;-0.015;0.010];

Dp.movingTheta0=deg2rad(0.20);
Dp.movingEOmega0=[0.010;-0.008;0.006];

Dp.initialAxis=[0.70;-0.20;0.40];
Dp.initialAxis=Dp.initialAxis/norm(Dp.initialAxis);

% Deterministic bounded torque residual for ISS/UUB testing [N m].
Dp.distAmp=[0.025;0.030;0.020];
Dp.distOmega=[0.60;0.40;0.80];
Dp.distPhase=[0.00;pi/2;0.30];

% Simulation settings.
Dp.Tconst=15.0;
Dp.Taggressive=120.0;
Dp.Twrap=50.0;
Dp.Tvarying=100.0;
Dp.tSteady=60.0;

Dp.RelTol=1e-9;
Dp.AbsTol=1e-11;
Dp.MaxStep=0.01;

% GNA damping-certificate scan.
Dp.scanNV=41;
Dp.scanNA=49;
end
