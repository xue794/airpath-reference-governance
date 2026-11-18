function M=final_parameter_manifest_rev1()
%FINAL_PARAMETER_MANIFEST_REV1 Canonical final simulation/controller manifest.
%
% This function deliberately reads the canonical module parameter functions.
% The full 6DOF integration must use these exact structs, without shadow
% retunes.

M.revision='Final Parameter Synchronization / Regression Freeze Rev.1';

M.G=smc_guidance_params();
M.E=wind_estimator_params();
M.R=rg_params();
M.A=viA_params();
M.B=viB_params();
M.F=viC_params();
M.Dp=viD_params();
M.P7=vii_params();

% Values that changed relative to the previously verified module bundle.
M.changed.viC_KomegaF_old=19.2*eye(3);
M.changed.viC_KomegaF_final=M.F.KomegaF;

% These actuator quantities are frozen numerical-study assumptions, not
% claimed T-2 hardware measurements.
M.studyAssumptions.deltaMin=M.P7.deltaMin;
M.studyAssumptions.deltaMax=M.P7.deltaMax;
M.studyAssumptions.tauDelta=M.P7.tauDelta;
M.studyAssumptions.deltaRateMax=M.P7.deltaRateMax;
M.studyAssumptions.tauT=M.P7.tauT;
M.studyAssumptions.TRateMax=M.P7.TRateMax;

% Legacy isolated VI-D disturbance-test level.  It is not a global
% physical full-loop hard bound.
M.dTauLegacyDiagnostic=M.P7.dTauBudget;
end
