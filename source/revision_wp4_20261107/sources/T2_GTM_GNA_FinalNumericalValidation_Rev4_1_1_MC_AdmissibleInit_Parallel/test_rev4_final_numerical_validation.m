function Results=test_rev4_final_numerical_validation()
%TEST_REV4_FINAL_NUMERICAL_VALIDATION One-command final validation acceptance.
%
% Order:
%   1. Rev.4 path/logic/parameter preflight.
%   2. Frozen Rev.1 module + complete-6DOF regression.
%   3. Rev.4 local macro-cascade sampled certificate.
%   4. Consolidated extended nonlinear performance metrics.
%   5. Export paper-facing ledgers and figures.

Pre=rev4_preflight();

outdir=fullfile(Pre.root,Pre.R.outputDir);
if ~exist(outdir,'dir'), mkdir(outdir); end

fprintf('\nREV.4 PHASE A -- FROZEN CONTROLLER / FULL-LOOP REGRESSION\n');
fprintf('=========================================================\n');
Reg=test_final_regression_freeze_rev1();
save(fullfile(outdir,'REV4_CHECKPOINT_PHASE_A.mat'),'Pre','Reg','-v7.3');

fprintf('\nREV.4 PHASE B -- LOCAL MACRO-CASCADE THEOREM INSTANCE\n');
fprintf('====================================================\n');
Cert=rev4_macro_certificate();
save(fullfile(outdir,'REV4_CHECKPOINT_PHASE_B.mat'),'Pre','Reg','Cert','-v7.3');

fprintf('\nREV.4 PHASE C -- EXTENDED NONLINEAR VALIDATION SUMMARY\n');
fprintf('=====================================================\n');
Ext=rev4_extended_metrics(Reg.CL,Cert);

assert(strcmp(Cert.status,'REV4_MACRO_CASCADE_CERTIFICATE_PASSED'), ...
    'Rev.4 local macro-cascade certificate did not pass.');

Results=struct();
Results.status='REV4_FINAL_NUMERICAL_VALIDATION_PASSED';
Results.Preflight=Pre;
Results.Regression=Reg;
Results.Certificate=Cert;
Results.Extended=Ext;

fprintf('\nREV.4 PHASE D -- EXPORT\n');
fprintf('======================\n');
Results.Files=rev4_export_results(Results);

fprintf('\nREV.4 FINAL NUMERICAL VALIDATION ACCEPTANCE\n');
fprintf('===========================================\n');
fprintf('All frozen module/full-loop regression assertions passed.\n');
fprintf('Local T-2 macro-cascade sampled certificate passed.\n');
fprintf('Extended varying/noisy/wrap/cold-start cases remain performance/robustness validations.\n');
fprintf('Canonical K_omega,f = %.1f I; no controller retuning was introduced.\n', ...
    Reg.CL.C.F.KomegaF(1,1));
fprintf('FINAL STATUS: %s\n',Results.status);
end
