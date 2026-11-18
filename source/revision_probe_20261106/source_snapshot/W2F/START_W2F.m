% START_W2F -- Paper 1 Aerospace publication computational timing benchmark.
% MATLAB R2018a. Deliberately SERIAL: do not start a parallel pool.
clear Report;
sourceRoot='C:\FWUAV\T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel';
Report=run_p1_aerospace_w2f(sourceRoot);
