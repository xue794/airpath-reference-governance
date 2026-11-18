% START_W2D -- Paper 1 -> Aerospace W2-D one-click targeted validation
% MATLAB R2018a. Keep this whole package together; do not copy isolated files.
packageRoot=fileparts(mfilename('fullpath'));
cd(packageRoot);
sourceRoot='C:\FWUAV\T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel';
Report=run_p1_aerospace_w2d(sourceRoot,'UseParallel',true,'NumWorkers',56); %#ok<NASGU>
