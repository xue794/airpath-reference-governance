function P7=vii_params()
%VII_PARAMS Section VII allocation + actuator-interface parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% IMPORTANT:
% - The +/-10 deg surface box is a conservative study box inside the
%   large-amplitude GNA validation excitation range.  It is NOT claimed
%   to be a documented GTM hardware stop.
% - Servo time constants/rate limits are module-test assumptions and remain
%   provisional until the final full-plant tuning.

R=rg_params();

P7.revision='Section VII Constrained Control Allocation + Actuator Interface Rev.0';

% Surface position box [aileron elevator rudder].
P7.deltaMin=deg2rad([-10;-10;-10]);
P7.deltaMax=deg2rad([ 10; 10; 10]);

% Surface servo model.
P7.tauDelta=[0.025;0.025;0.025];       % s
P7.deltaRateMax=deg2rad([120;120;120]); % rad/s

% Propulsor interface.
P7.Tmin=R.Tmin;
P7.Tmax=R.Tmax;
P7.tauT=0.08;                         % s
P7.TRateMax=100.0;                    % N/s

% Weighted box-constrained least-squares allocator.
P7.Wtau=eye(3);
P7.allocTol=1e-11;

% VI-D disturbance-budget interface.
P7.dTauBudget=0.0405;                 % N m, legacy VI-D diagnostic level only
P7.tSteady=60.0;

% Stress request for saturation / active-set unit test.
P7.stressTau=[20;8;15];               % N m

% Simulation settings.
P7.Taggressive=120.0;
P7.Tvarying=100.0;
P7.Twrap=50.0;
P7.RelTol=1e-9;
P7.AbsTol=1e-11;
P7.MaxStep=0.005;

% Deterministic effectiveness scan.
P7.scanNV=41;
P7.scanNA=49;
end
