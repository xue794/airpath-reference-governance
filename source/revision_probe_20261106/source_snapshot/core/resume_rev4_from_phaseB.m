function Results4=resume_rev4_from_phaseB()
%RESUME_REV4_FROM_PHASEB Resume Rev.4 after a Phase-C/D reporting hotfix.
%
% Requires:
%   rev4_outputs/REV4_CHECKPOINT_PHASE_B.mat
%
% This avoids repeating the frozen regression and macro-certificate solve
% when those phases have already completed successfully under Rev.4.0.1.

root=fileparts(mfilename('fullpath'));
R=rev4_params();
ck=fullfile(root,R.outputDir,'REV4_CHECKPOINT_PHASE_B.mat');

if ~exist(ck,'file')
    error('Rev4:Checkpoint','Phase-B checkpoint not found: %s',ck);
end

L=load(ck,'Pre','Reg','Cert');
Pre=L.Pre;
Reg=L.Reg;
Cert=L.Cert;

fprintf('\nREV.4 RESUME FROM PHASE B CHECKPOINT\n');
fprintf('===================================\n');
fprintf('Loaded successful frozen regression and macro-certificate checkpoint.\n');

Ext=rev4_extended_metrics(Reg.CL,Cert);

assert(strcmp(Cert.status,'REV4_MACRO_CASCADE_CERTIFICATE_PASSED'), ...
    'Checkpoint macro-cascade certificate status is not PASS.');

Results4=struct();
Results4.status='REV4_FINAL_NUMERICAL_VALIDATION_PASSED';
Results4.Preflight=Pre;
Results4.Regression=Reg;
Results4.Certificate=Cert;
Results4.Extended=Ext;
Results4.Files=rev4_export_results(Results4);

fprintf('\nREV.4 FINAL NUMERICAL VALIDATION ACCEPTANCE\n');
fprintf('===========================================\n');
fprintf('Resumed from successful Phase-B checkpoint.\n');
fprintf('FINAL STATUS: %s\n',Results4.status);

assignin('base','Results4',Results4);
end
