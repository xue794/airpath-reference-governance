function K=w2d_case_info(mode,C,O)
%W2D_CASE_INFO Predeclared targeted deterministic cases.
if nargin<3||isempty(O),O=w2d_options();end
K=struct();K.mode=mode;K.windMode='still';K.addNoise=false;K.coldWindStart=false;
K.T=O.CaseDuration_s;K.tSteady=O.LateStart_s;K.RampStart_s=O.RampStart_s;K.RampDuration_s=O.RampDuration_s;K.PulseHold_s=O.PulseHold_s;
switch lower(mode)
 case 'target_tmin'
  K.TargetAirPath=O.TargetTmin;K.ExpectedRowIndex=16;K.ExpectedRow='T_min';
 case 'target_tmax'
  K.TargetAirPath=O.TargetTmax;K.ExpectedRowIndex=17;K.ExpectedRow='T_max';
 case 'target_nplus'
  K.TargetAirPath=O.TargetNplus;K.ExpectedRowIndex=15;K.ExpectedRow='N_plus';
 otherwise
  error('W2D:Mode','Unknown W2-D mode: %s',mode);
end
K.InitialAirPath=[C.R.Vb;C.R.gammaB;0];
end
