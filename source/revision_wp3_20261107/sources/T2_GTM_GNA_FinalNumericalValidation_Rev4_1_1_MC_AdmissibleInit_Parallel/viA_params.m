function A = viA_params()
%VIA_PARAMS Section VI-A air-path inversion / ideal realization parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% These values are retained as the FINAL synchronized Rev.1 air-path gains.

A.revision = 'Section VI-A Air-Path Force Inversion and Ideal Realization Rev.0';

A.kV = 0.80;               % 1/s
A.kGamma = 1.20;           % 1/s
A.kChi = 1.00;             % 1/s

% Initial path tracking errors [m/s; rad; rad].
A.e0 = [0.30; deg2rad(0.30); deg2rad(-0.50)];
A.e0Wrap = [0.20; deg2rad(0.20); deg2rad(-0.40)];

A.Tstill = 60;
A.Tconstant = 60;
A.Tvarying = 100;
A.Taggressive = 120;
A.Twrap = 50;
A.sample_dt = 0.05;

A.RelTol = 1e-9;
A.AbsTol = 1e-11;

% Steady window for the time-varying-wind UUB check.
A.tUUB0 = 60;
end
