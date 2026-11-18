function P=rev4_preflight()
%REV4_PREFLIGHT Path, toolbox, parameter, and logic-lock audit.

root=fileparts(mfilename('fullpath'));
R=rev4_params();
C=cl6_params();
M=final_parameter_manifest_rev1();

fprintf('\nREV.4 FINAL NUMERICAL VALIDATION PREFLIGHT\n');
fprintf('==========================================\n');
fprintf('root = %s\n',root);

mustBeLocal={ ...
    'rev4_params', ...
    'rev4_macro_certificate', ...
    'rev4_extended_metrics', ...
    'rev4_export_results', ...
    'test_rev4_final_numerical_validation', ...
    'test_final_regression_freeze_rev1', ...
    'simulate_complete_6dof_case', ...
    'viii_reduced_vector_field', ...
    'viii_lift_reduced', ...
    'viii_reduced_from_full'};

for k=1:numel(mustBeLocal)
    p=which(mustBeLocal{k});
    if isempty(p)
        error('Rev4:Path','Required function is not on the path: %s',mustBeLocal{k});
    end
    if length(p)<length(root) || ~strcmpi(p(1:length(root)),root)
        error('Rev4:Shadow','Function is shadowed outside the Rev.4 root: %s -> %s',mustBeLocal{k},p);
    end
end

hasLMI=(exist('setlmis','file')==2) && (exist('lmivar','file')==2) && ...
       (exist('lmiterm','file')==2) && (exist('feasp','file')==2) && ...
       (exist('dec2mat','file')==2);

fprintf('PRE-P1 function shadow audit = PASS\n');
fprintf('PRE-P2 Robust-Control LMI functions available = %d\n',hasLMI);

if ~hasLMI
    error('Rev4:LMI','Rev.4 final macro certificate requires the Robust Control Toolbox LMI functions.');
end

resF=norm(C.F.KomegaF-M.F.KomegaF,'fro');
assert(abs(M.F.kRf-256)<1e-14 && ...
       norm(M.F.KomegaF-7*eye(3),'fro')<1e-14 && ...
       abs(M.F.cF-0.5)<1e-14 && resF<1e-14,'Rev4 parameter lock failed.');

fprintf('PRE-P3 canonical K_omega,f = %.1f I, sync residual = %.3e\n', ...
    M.F.KomegaF(1,1),resF);
fprintf('PRE-P4 certified-small validation scale = %.6f (validation only)\n',R.certScale);
fprintf('PRE-P5 theory reference = %s\n',R.theoryFile);
fprintf('PRE-P6 theory SHA-256 = %s\n',R.theorySHA256);
fprintf('PRE-P7 logic lock: final certificate = flight-control macro + autonomous wind-estimator cascade\n');
fprintf('PRE-P8 old eight-block scalar small-gain tests are historical diagnostics only\n');

P=struct('root',root,'R',R,'C',C,'manifest',M,'hasLMI',hasLMI);
end
