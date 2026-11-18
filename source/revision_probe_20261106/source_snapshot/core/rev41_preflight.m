function P=rev41_preflight(R)
%REV41_PREFLIGHT Rev.4.1 path, logic, and frozen-controller audit.
if nargin<1 || isempty(R), R=rev41_mc_params(); end

root=fileparts(mfilename('fullpath'));
C=cl6_params();
M=final_parameter_manifest_rev1();

fprintf('\nREV.4.1.1 MONTE-CARLO PREFLIGHT\n');
fprintf('==============================\n');
fprintf('root = %s\n',root);

must={ ...
    'run_rev41_monte_carlo', ...
    'rev41_generate_samples', ...
    'rev41_apply_sample', ...
    'rev41_run_one', ...
    'simulate_complete_6dof_case', ...
    'cl6_case_info', ...
    'cl6_wind_profile', ...
    'wind_estimator_test_noise'};

for k=1:numel(must)
    q=which(must{k});
    if isempty(q)
        error('Rev41:Path','Required function missing: %s',must{k});
    end
    if length(q)<length(root) || ~strcmpi(q(1:length(root)),root)
        error('Rev41:Shadow','Function shadowed outside Rev.4.1 root: %s -> %s',must{k},q);
    end
end

resF=norm(C.F.KomegaF-M.F.KomegaF,'fro');
assert(abs(M.F.kRf-256)<1e-14 && ...
       norm(M.F.KomegaF-7*eye(3),'fro')<1e-14 && ...
       abs(M.F.cF-0.5)<1e-14 && resF<1e-14, ...
       'Rev.4.1 frozen-controller lock failed.');

fprintf('MC-PRE1 function shadow audit = PASS\n');
fprintf('MC-PRE2 canonical K_omega,f = %.1f I, sync residual = %.3e\n', ...
    M.F.KomegaF(1,1),resF);
fprintf('MC-PRE3 trials / seed = %d / %d\n',R.N,R.seed);
fprintf('MC-PRE4 study-defined uncertainty set; no controller retuning\n');
fprintf('MC-PRE5 Monte Carlo is extended empirical robustness, not theorem proof\n');
fprintf('MC-PRE6 initialization strategy = %s\n',R.initStrategy);
fprintf('MC-PRE7 all uncertainty draws are retained; no rejection sampling\n');

P=struct('root',root,'R',R,'Cnom',C,'manifest',M);
end
