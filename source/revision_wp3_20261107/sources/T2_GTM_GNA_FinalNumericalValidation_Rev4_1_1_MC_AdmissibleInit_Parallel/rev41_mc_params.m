function R=rev41_mc_params()
%REV41_MC_PARAMS Frozen randomized robustness-study definition.
%
% IMPORTANT:
% These distributions are study-defined robustness distributions. They are
% NOT identified probability distributions for the physical GTM vehicle.
% Consequently, empirical failure rates apply only to this stated sampling
% model and are not claimed as real-world reliability probabilities.

R.revision='Rev.4.1.1 Monte-Carlo Robustness Assessment -- Admissible Initialization';
R.N=500;
R.seed=20260819;
R.mode='monte_carlo';
R.T=120.0;
R.tSteady=60.0;
R.checkpointEvery=10;
R.resume=true;

% Reference-governor initialization.
% All randomized trials start at the strict backup-feasible filter state;
% the random raw command itself is NOT altered.
R.initStrategy='backup_feasible';

% Controller is frozen. Only plant/environment/actuator/initial conditions vary.
R.controllerFrozen=true;

% Rigid-body plant uncertainty: independent uniform multiplicative scales.
R.massFrac=0.05;       % mass in [0.95,1.05] nominal
R.inertiaFrac=0.05;    % Ixx,Iyy,Izz,Ixz independently +/-5%

% Aerodynamic model mismatch: every GNA coefficient gets an independent
% multiplicative factor in [0.95,1.05].  This is an epistemic stress ensemble,
% not a parameter-identification covariance model.
R.aeroFrac=0.05;

% Actuator first-order time constants only; hard boxes and rate limits stay fixed.
R.servoTauFrac=0.20;   % tau_delta independently +/-20%
R.thrustTauFrac=0.20;  % tau_T +/-20%

% Wind family around the deterministic varying-wind study case.
R.windMean0=[6;-3;0];
R.windMeanHalfRange=[2;1.5;0.5];
R.windAmp0=[2;1.5;0.8];
R.windAmpScaleRange=[0.75 1.25];
R.windOmega0=[0.25;0.40;0.20];
R.windOmegaScaleRange=[0.80 1.20];

% Sensor-error waveform. Scale <=1 preserves all original Rev.4 error bounds.
R.noiseScaleRange=[0.50 1.00];
R.noiseFreqScaleRange=[0.80 1.20];

% Initial-condition perturbations.
R.initPosHalfRange=[2;2;1];                 % m
R.initVgBHalfRange=[0.30;0.30;0.20];       % m/s
R.initAttitudeMax=deg2rad(0.20);            % principal angle
R.initOmegaHalfRange=deg2rad([0.50;0.50;0.50]); % rad/s
R.initWhatHalfRange=[0.50;0.50;0.30];       % m/s

% Numerical/constraint acceptance tolerances.
R.phiTol=1e-7;
R.boxTol=1e-7;
R.rateTol=1e-7;
R.qNormTol=1e-7;
R.stateTol=1e-8;

% Certified-domain diagnostics (reported separately from hard feasibility).
R.totalAttitudeTube=deg2rad(0.70);
R.betaTube=deg2rad(0.70);
R.kappaCert=0.985622042;

% Export.
R.outputDir='rev411_outputs';
R.makePlots=true;
end
