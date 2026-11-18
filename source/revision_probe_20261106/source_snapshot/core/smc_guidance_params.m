function G = smc_guidance_params()
%SMC_GUIDANCE_PARAMS Section IV integral-error SMC guidance Rev.0 parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% These values are retained as the FINAL synchronized Rev.1 simulation/controller gains.

G.revision = 'Section IV Wind-Compensated Integral-Error SMC Guidance Rev.0';

% s = e_r + Lambda z_e
G.Lambda = diag([0.40 0.50 0.60]);       % [1/s]
G.Ks     = diag([0.80 0.90 1.00]);       % [1/s]
G.Kc     = [1.20; 1.20; 0.80];           % [m/s]
G.eps    = [0.50; 0.50; 0.40];           % [m]

% Initial position error used in module tests.
G.e0 = [30; -20; 10];                    % [m]
G.ze0 = [0;0;0];                         % [m*s]

% Test durations and steady-state windows.
G.t_still = 60.0;
G.t_constant = 60.0;
G.t_varying = 120.0;
G.t_mismatch = 80.0;
G.t_bound = 140.0;
G.fit_window = [60 120];
G.bound_window = [80 140];

% ODE tolerances chosen to expose algebraic identities numerically.
G.RelTol = 1e-10;
G.AbsTol = 1e-12;
end
