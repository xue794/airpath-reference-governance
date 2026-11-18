% START_W2C  Paper 1 -> Aerospace, W2-C one-click deterministic comparison.
% MATLAB R2018a. Run this script from the extracted W2-C package folder.
clear Report;
pkgRoot=fileparts(mfilename('fullpath'));
cd(pkgRoot);
sourceRoot='C:\FWUAV\T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel';
Report=run_p1_aerospace_w2c(sourceRoot,'UseParallel',true,'NumWorkers',56);
