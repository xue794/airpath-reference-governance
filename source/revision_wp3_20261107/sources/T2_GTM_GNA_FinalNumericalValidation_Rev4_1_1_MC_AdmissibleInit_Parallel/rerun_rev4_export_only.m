function Files=rerun_rev4_export_only()
%RERUN_REV4_EXPORT_ONLY Recreate Rev.4 figures/ledgers after exporter-only hotfix.
%
% This function does not rerun the controller simulations or the macro LMI.
% It reuses the already accepted Results4 structure, either from the base
% workspace or from rev4_outputs/Rev4Results.mat.

root=fileparts(mfilename('fullpath'));
R=rev4_params();
matFile=fullfile(root,R.outputDir,'Rev4Results.mat');

hasBase=evalin('base','exist(''Results4'',''var'')==1');

if hasBase
    Results4=evalin('base','Results4');
else
    if ~exist(matFile,'file')
        error('Rev4:ExportOnly','Neither base-workspace Results4 nor %s exists.',matFile);
    end
    L=load(matFile,'Results');
    Results4=L.Results;
end

if ~isfield(Results4,'status') || ...
        ~strcmp(Results4.status,'REV4_FINAL_NUMERICAL_VALIDATION_PASSED')
    error('Rev4:ExportOnly','Stored Results4 is not an accepted Rev.4 result.');
end

fprintf('\nREV.4.0.2 EXPORT-ONLY HOTFIX\n');
fprintf('============================\n');
fprintf('Reusing accepted Rev.4 numerical results; simulations/LMI are not rerun.\n');

Files=rev4_export_results(Results4);
Results4.Files=Files;
assignin('base','Results4',Results4);

fprintf('EXPORT-ONLY STATUS: REV4_0_2_EXPORT_PASSED\n');
end
