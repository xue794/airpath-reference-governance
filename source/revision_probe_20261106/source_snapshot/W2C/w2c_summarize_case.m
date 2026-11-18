function E=w2c_summarize_case(S,O)
%W2C_SUMMARIZE_CASE Flat scalar summary for CSV export and gating.
Y=S.Y;erNorm=sqrt(sum(Y.er.^2,2));late=S.t>=S.Case.tSteady;if ~any(late),late=true(size(S.t));end
n=numel(S.t);VB=zeros(n,1);for k=1:n,xi=rg_backup_state(Y.xf(k,:).',S.C.R);VB(k)=xi.'*S.C.R.PB*xi;end
windErr=sqrt(sum(Y.tildew.^2,2));momRes=sqrt(sum(Y.dTauPhysical.^2,2));
E=struct();E.Method=S.Method;E.Mode=S.mode;E.Status='COMPLETED';E.ErrorId='';E.ErrorMessage='';E.rhoB=S.C.R.rhoB;
E.IntegrationWall_s=S.IntegrationWall_s;E.OutputWall_s=S.OutputWall_s;
E.PosRMS_Full_m=sqrt(mean(erNorm.^2));E.PosPeak_Full_m=max(erNorm);E.PosRMS_Late_m=sqrt(mean(erNorm(late).^2));E.PosPeak_Late_m=max(erNorm(late));
E.GovInterventionPct=100*mean(Y.intervened);E.MeanCorrectionNorm=mean(Y.correctionNorm);E.MedianCorrectionNorm=median(Y.correctionNorm);E.MaxCorrectionNorm=max(Y.correctionNorm);
E.MeanObjective=mean(Y.objective);E.MedianObjective=median(Y.objective);E.MaxObjective=max(Y.objective);
if any(Y.intervened),E.MeanObjectiveActive=mean(Y.objective(Y.intervened));E.MedianObjectiveActive=median(Y.objective(Y.intervened));else,E.MeanObjectiveActive=0;E.MedianObjectiveActive=0;end
if strcmpi(S.Method,'Proposed'),E.MinScalarLambda=min(Y.scalarLambda);else,E.MinScalarLambda=NaN;end
E.MaxVB=max(VB);E.MaxVB_over_rhoB=max(VB)/S.C.R.rhoB;E.MinBackupMargin=min(S.C.R.rhoB-VB);
E.MinAppliedPhi=min(Y.phi(:));E.MaxNormalizedPhiViolation=max(Y.normalizedPhiViolation);
E.MaxAttitudeTube_deg=rad2deg(max(Y.thetaF+Y.thetaR));E.LateMaxWindEstErr_mps=max(windErr(late));
E.ActualVmin_mps=min(Y.V);E.ActualVmax_mps=max(Y.V);E.ActualMaxAbsGamma_deg=rad2deg(max(abs(Y.gamma)));
E.ActualAlphaMin_deg=rad2deg(min(Y.alpha));E.ActualAlphaMax_deg=rad2deg(max(Y.alpha));E.ActualMaxAbsBeta_deg=rad2deg(max(abs(Y.beta)));E.MinKappa=min(Y.kappa);
E.ActualTmin_N=min(Y.Tact);E.ActualTmax_N=max(Y.Tact);E.MaxAbsSurface_deg=max(max(abs(rad2deg(Y.deltaAct))));E.MaxAbsSurfaceRate_degps=max(max(abs(rad2deg(Y.deltaDot))));
E.AllocatorActiveSamples=sum(any(Y.active>0,2));E.MaxMomentResidual_Nm=max(momRes);E.MaxThrustTrackingErr_N=max(abs(Y.Tact-Y.Tcmd));E.MaxThrustRate_Nps=max(abs(Y.TDot));
E.MaxSlidingIdentityResidual=max(Y.slidingIdentityResidual);E.MaxWindKinematicResidual=max(Y.windKinematicResidual);
E.QPEnumerationFallbackSamples=sum(Y.qpEnumerationFallback);E.QPFeasiblePathSamples=sum(strcmp(Y.qpPath,'feasible'));E.QPSingleHalfspaceSamples=sum(strcmp(Y.qpPath,'single-halfspace'));
E.QPPrimalActiveSetSamples=sum(strcmp(Y.qpPath,'primal-active-set'));E.QPEnumerationSamples=sum(strncmp(Y.qpPath,'enumeration',11));
if strcmpi(S.Method,'MatchedQP')
 E.QPMaxPrimal=max(Y.qpPrimal);E.QPMaxStationarity=max(Y.qpStationarity);E.QPMaxDualViolation=max(Y.qpDualViolation);E.QPMaxComplementarity=max(Y.qpComplementarity);
else
 E.QPMaxPrimal=NaN;E.QPMaxStationarity=NaN;E.QPMaxDualViolation=NaN;E.QPMaxComplementarity=NaN;
end
% Physical box audit (model-defined study limits, not hardware certification).
E.ActuatorBoxPass=E.ActualTmin_N>=S.C.P7.Tmin-1e-9 && E.ActualTmax_N<=S.C.P7.Tmax+1e-9 && ...
 E.MaxAbsSurface_deg<=rad2deg(max(S.C.P7.deltaMax))+1e-7 && E.MaxAbsSurfaceRate_degps<=rad2deg(max(S.C.P7.deltaRateMax))+1e-7 && E.MaxThrustRate_Nps<=S.C.P7.TRateMax+1e-7;
E.ReferenceGatePass=E.MaxNormalizedPhiViolation<=O.NormalizedPhiTol && E.MaxVB_over_rhoB<=1+O.BackupTol;
if strcmpi(S.Method,'MatchedQP')
 E.QPKKTGatePass=E.QPMaxPrimal<=O.KKTGateTol && E.QPMaxStationarity<=O.KKTGateTol && E.QPMaxDualViolation<=O.KKTGateTol && E.QPMaxComplementarity<=O.KKTGateTol;
else,E.QPKKTGatePass=true;end
end
