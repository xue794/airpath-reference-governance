%% run_transparent_governor_baseline_R2018a
% Controlled ablation for the CJA manuscript.
%
% Purpose:
%   Generate a fair "transparent governor" baseline for the accepted
%   time-varying-wind case. Everything is held fixed except the certified
%   scalar interpolation:
%
%       proposed:    u* = u_B + lambda*(u0-u_B),  0 <= lambda* <= 1
%       baseline:    u* = u0,                     lambda* = 1
%
% The same:
%   - T-2/GNA plant
%   - wind realization
%   - initial condition
%   - wind estimator
%   - SMC guidance
%   - nominal second-order reference generator
%   - air-path realization
%   - PHS rotational controller
%   - allocator
%   - physical actuator dynamics
%   - controller gains
% are retained. No retuning is performed.
%
% IMPORTANT:
%   Place this file and the other ablation helper .m files in the SAME
%   directory as the frozen final numerical-validation code, or add both
%   directories to the MATLAB path.
%
% MATLAB R2018a compatible.

clc; close all;

MODE = 'varying';
OUT_FILE = 'TransparentGovernorBaseline_varying.mat';

required = { ...
    'cl6_params', ...
    'cl6_case_info', ...
    'cl6_initial_state', ...
    't2_gna_6dof_quat', ...
    'cl6_plant_input', ...
    'rg_raw_command_chart', ...
    'rg_nominal_control', ...
    'rg_constraint_affine', ...
    'rg_backup_control'};

for k = 1:numel(required)
    if exist(required{k},'file') ~= 2
        error(['Required frozen-code function not found: %s\n' ...
               'Run this script from the final numerical-validation code directory, ' ...
               'or add that directory to the MATLAB path.'], required{k});
    end
end

C = cl6_params();

fprintf('\n============================================================\n');
fprintf('Transparent-governor controlled ablation\n');
fprintf('Mode: %s\n', MODE);
fprintf('Only certified scalar interpolation is disabled (lambda*=1).\n');
fprintf('No controller/gain retuning is performed.\n');
fprintf('============================================================\n');

try
    Baseline = simulate_complete_6dof_case_transparent_R2018a(MODE, C);
catch ME
    fprintf(2,'\nTransparent-governor simulation did not complete.\n');
    fprintf(2,'Identifier: %s\n', ME.identifier);
    fprintf(2,'Message   : %s\n\n', ME.message);
    fprintf(2,['This outcome may itself indicate that the transparent reference leaves ' ...
               'the implementation domain. Do NOT weaken or clip the baseline silently.\n']);
    rethrow(ME);
end

Baseline.ablation = struct();
Baseline.ablation.name = 'transparent_governor';
Baseline.ablation.description = ...
    'Same frozen complete-loop architecture and gains; certified scalar governor interpolation disabled by enforcing lambda*=1 and u*=u0.';
Baseline.ablation.certificationEnabled = false;
Baseline.ablation.lambdaForced = 1;
Baseline.ablation.retuned = false;
Baseline.ablation.mode = MODE;

save(OUT_FILE, 'Baseline', '-v7');

fprintf('\nBaseline generated successfully.\n');
fprintf('Saved to: %s\n', fullfile(pwd,OUT_FILE));
fprintf('Samples: %d\n', numel(Baseline.t));
fprintf('min feasibility margin: %.12g\n', min(Baseline.Y.minPhi));
fprintf('max position-error norm: %.6f m\n', max(sqrt(sum(Baseline.Y.er.^2,2))));
fprintf('RMS position-error norm: %.6f m\n\n', sqrt(mean(sum(Baseline.Y.er.^2,2))));
