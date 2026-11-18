function E=w2e_summarize_run(S,Samp,Rmc,O)
%W2E_SUMMARIZE_RUN Complete per-run scalar/vector diagnostics and frozen gates.
Y=S.Y;C=S.C;E=w2e_empty_run();E.Index=Samp.index;E.Method=S.Method;E.Status='COMPLETED';E.Completed=true;E.FailureClass='NONE';E.rhoB=C.R.rhoB;E.InitialVB_over_rhoB=S.I0.VB0_over_rho;E.IntegrationWall_s=S.IntegrationWall_s;E.OutputWall_s=S.OutputWall_s;
erNorm=sqrt(sum(Y.er.^2,2));late=S.t>=S.Case.tSteady;if ~any(late),late=true(size(S.t));end
E.PosRMSFull_m=sqrt(mean(erNorm.^2));E.PosPeakFull_m=max(erNorm);E.PosRMSLate_m=sqrt(mean(erNorm(late).^2));E.PosPeakLate_m=max(erNorm(late));E.LateErLinf_m=max(abs(Y.er(late,:)),[],1);
windErr=sqrt(sum(Y.tildew.^2,2));momRes=sqrt(sum(Y.dTauPhysical.^2,2));surfForce=sqrt(sum(Y.fSurface.^2,2));E.LateWindErrorMax_mps=max(windErr(late));
E.GovInterventionPct=100*mean(Y.intervened);E.InterventionSamples=sum(Y.intervened);E.MeanCorrectionNorm=mean(Y.correctionNorm);E.MedianCorrectionNorm=median(Y.correctionNorm);E.MaxCorrectionNorm=max(Y.correctionNorm);
E.MeanObjective=mean(Y.objective);E.MedianObjective=median(Y.objective);E.MaxObjective=max(Y.objective);E.IntegratedObjective_s=trapz(S.t,Y.objective);
if any(Y.intervened),E.MeanObjectiveActive=mean(Y.objective(Y.intervened));else,E.MeanObjectiveActive=0;end
if strcmpi(S.Method,'Proposed'),E.MinScalarLambda=min(Y.scalarLambda);end
n=numel(S.t);VB=zeros(n,1);for k=1:n,xi=rg_backup_state(Y.xf(k,:).',C.R);VB(k)=xi.'*C.R.PB*xi;end
E.MaxVB_over_rhoB=max(VB)/C.R.rhoB;E.MinBackupMargin=min(C.R.rhoB-VB);E.MinAppliedPhi=min(Y.phi(:));E.MaxNormalizedPhiViolation=max(Y.normalizedPhiViolation);
E.MaxAttitudeTube_deg=rad2deg(max(Y.thetaF+Y.thetaR));E.ActualVmin_mps=min(Y.V);E.ActualVmax_mps=max(Y.V);E.ActualMaxAbsGamma_deg=rad2deg(max(abs(Y.gamma)));E.ActualAlphaMin_deg=rad2deg(min(Y.alpha));E.ActualAlphaMax_deg=rad2deg(max(Y.alpha));E.ActualMaxAbsBeta_deg=rad2deg(max(abs(Y.beta)));E.MinKappa=min(Y.kappa);
E.ActualTmin_N=min(Y.Tact);E.ActualTmax_N=max(Y.Tact);E.MaxThrustRate_Nps=max(abs(Y.TDot));dm=max(abs(rad2deg(Y.deltaAct)),[],1);dr=max(abs(rad2deg(Y.deltaDot)),[],1);E.Delta1Max_deg=dm(1);E.Delta2Max_deg=dm(2);E.Delta3Max_deg=dm(3);E.DeltaRate1Max_degps=dr(1);E.DeltaRate2Max_degps=dr(2);E.DeltaRate3Max_degps=dr(3);
E.AllocatorActiveSamples=sum(any(Y.active>0,2));E.MaxMomentResidual_Nm=max(momRes);E.LateMomentResidual_Nm=max(momRes(late));E.MaxSurfaceForce_N=max(surfForce);E.MaxSlidingIdentityResidual=max(Y.slidingIdentityResidual);E.MaxWindKinematicResidual=max(Y.windKinematicResidual);E.MaxQNormError=max(abs(Y.qNorm-1));E.MaxQrNormError=max(abs(Y.qrNorm-1));
% Per-row records. "Role" means scalar limiter for Proposed and positive-dual support for QP.
if strcmpi(S.Method,'Proposed'),E.RoleCount=sum(Y.scalarLimiter,1);else,E.RoleCount=sum(Y.qpSupport,1);end
E.BindingCount=sum(Y.binding,1);E.AffineMin=min(Y.phi,[],1);E.NominalAffineMin=min(Y.phi0,[],1);E.StateMarginMin=min(Y.stateMargins,[],1);
if strcmpi(S.Method,'Proposed')
 Rsel=Y.scalarLimiter;
else
 Rsel=Y.qpSupport;
end
E.NonbackupPhysicalRoleSamples=sum(any(Rsel(:,11:17),2));
E.ThrustRoleSamples=sum(any(Rsel(:,16:17),2));
E.TminRoleSamples=E.RoleCount(16);E.TmaxRoleSamples=E.RoleCount(17);E.BackupRoleSamples=E.RoleCount(18);
E.QPEnumerationFallbackSamples=sum(Y.qpEnumerationFallback);E.QPFeasiblePathSamples=sum(strcmp(Y.qpPath,'feasible'));E.QPSingleHalfspaceSamples=sum(strcmp(Y.qpPath,'single-halfspace'));E.QPPrimalActiveSetSamples=sum(strcmp(Y.qpPath,'primal-active-set'));E.QPEnumerationSamples=sum(strncmp(Y.qpPath,'enumeration',11));
if strcmpi(S.Method,'MatchedQP')
 E.QPMaxPrimal=max(Y.qpPrimal);E.QPMaxStationarity=max(Y.qpStationarity);E.QPMaxDualViolation=max(Y.qpDualViolation);E.QPMaxComplementarity=max(Y.qpComplementarity);
 E.QPKKTGatePass=E.QPMaxPrimal<=O.KKTGateTol && E.QPMaxStationarity<=O.KKTGateTol && E.QPMaxDualViolation<=O.KKTGateTol && E.QPMaxComplementarity<=O.KKTGateTol;
else
 E.QPKKTGatePass=true;
end
% Frozen Rev.4.1.1 hard-feasibility criteria, plus separately reported backup/normalized/KKT gates.
E.ReferenceAffinePass=E.MinAppliedPhi>=-Rmc.phiTol;E.ReferenceNormalizedPass=E.MaxNormalizedPhiViolation<=O.NormalizedPhiTol;E.BackupSetPass=E.MaxVB_over_rhoB<=1+O.BackupTol;
E.AirspeedPass=E.ActualVmin_mps>=C.R.Vmin-Rmc.stateTol && E.ActualVmax_mps<=C.R.Vmax+Rmc.stateTol;E.GammaPass=E.ActualMaxAbsGamma_deg<=rad2deg(C.R.gammaMax)+1e-8;E.AlphaPass=E.ActualAlphaMin_deg>=rad2deg(C.R.alphaPhysicalMin)-1e-8 && E.ActualAlphaMax_deg<=rad2deg(C.R.alphaPhysicalMax)+1e-8;
dmin=rad2deg(C.P7.deltaMin(:)).';dmax=rad2deg(C.P7.deltaMax(:)).';rmax=rad2deg(C.P7.deltaRateMax(:)).';
deltaDeg=rad2deg(Y.deltaAct);
E.SurfaceBoxPass=all(max(deltaDeg,[],1)<=dmax+rad2deg(Rmc.boxTol)) && all(min(deltaDeg,[],1)>=dmin-rad2deg(Rmc.boxTol));
E.SurfaceRatePass=all([E.DeltaRate1Max_degps E.DeltaRate2Max_degps E.DeltaRate3Max_degps]<=rmax+rad2deg(Rmc.rateTol));E.ThrustBoxPass=E.ActualTmin_N>=C.P7.Tmin-Rmc.boxTol && E.ActualTmax_N<=C.P7.Tmax+Rmc.boxTol;E.ThrustRatePass=E.MaxThrustRate_Nps<=C.P7.TRateMax+Rmc.rateTol;E.QuaternionPass=E.MaxQNormError<=Rmc.qNormTol && E.MaxQrNormError<=Rmc.qNormTol;E.KappaPositivePass=E.MinKappa>0;
E.ActuatorBoxPass=E.SurfaceBoxPass&&E.SurfaceRatePass&&E.ThrustBoxPass&&E.ThrustRatePass;
E.HardFeasible=E.Completed&&E.ReferenceAffinePass&&E.AirspeedPass&&E.GammaPass&&E.AlphaPass&&E.SurfaceBoxPass&&E.SurfaceRatePass&&E.ThrustBoxPass&&E.ThrustRatePass&&E.QuaternionPass&&E.KappaPositivePass;
E.CorePass=E.HardFeasible&&E.ReferenceNormalizedPass&&E.BackupSetPass&&E.QPKKTGatePass;
E.AttitudeTubePass=deg2rad(E.MaxAttitudeTube_deg)<=Rmc.totalAttitudeTube+Rmc.stateTol;E.BetaTubePass=deg2rad(E.ActualMaxAbsBeta_deg)<=Rmc.betaTube+Rmc.stateTol;E.KappaCertPass=E.MinKappa>=Rmc.kappaCert-Rmc.stateTol;E.CertifiedDomain=E.CorePass&&E.AttitudeTubePass&&E.BetaTubePass&&E.KappaCertPass;
end
