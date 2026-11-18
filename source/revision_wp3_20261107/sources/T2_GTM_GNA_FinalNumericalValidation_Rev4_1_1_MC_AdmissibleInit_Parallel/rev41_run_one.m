function Q=rev41_run_one(Samp,R,C0)
%REV41_RUN_ONE Execute and summarize one randomized complete-6DOF trial.

if nargin<3 || isempty(C0), C0=cl6_params(); end

Q=local_empty_run();
Q.index=Samp.index;

try
    C=rev41_apply_sample(C0,Samp,R);
    Sim=simulate_complete_6dof_case(R.mode,C);

    if any(~isfinite(Sim.x(:)))
        error('Rev41:Nonfinite','Nonfinite state encountered.');
    end

    Y=Sim.Y;
    late=Sim.t>=Sim.Case.tSteady;
    if ~any(late), late=true(size(Sim.t)); end

    er=Y.er;
    posNorm=sqrt(sum(er.^2,2));
    windErr=sqrt(sum(Y.tildew.^2,2));
    dTau=sqrt(sum(Y.dTauPhysical.^2,2));
    dF=sqrt(sum(Y.fSurface.^2,2));

    Q.solverOK=true;
    Q.lateErLinf=max(abs(er(late,:)),[],1);
    Q.latePosNormMax=max(posNorm(late));
    Q.latePosNormRMS=sqrt(mean(posNorm(late).^2));
    Q.maxPosNorm=max(posNorm);
    Q.lateWindErrorMax=max(windErr(late));

    Q.minLambda=min(Y.lambda);
    Q.minPhi=min(Y.minPhi);
    Q.maxTubeDeg=rad2deg(max(Y.thetaF+Y.thetaR));

    Q.Vmin=min(Y.V);
    Q.Vmax=max(Y.V);
    Q.gammaAbsMaxDeg=rad2deg(max(abs(Y.gamma)));
    Q.alphaMinDeg=rad2deg(min(Y.alpha));
    Q.alphaMaxDeg=rad2deg(max(Y.alpha));
    Q.betaAbsMaxDeg=rad2deg(max(abs(Y.beta)));
    Q.kappaMin=min(Y.kappa);

    Q.Tmin=min(Y.Tact);
    Q.Tmax=max(Y.Tact);
    Q.maxThrustRate=max(abs(Y.TDot));
    Q.maxDeltaActDeg=max(abs(rad2deg(Y.deltaAct)),[],1);
    Q.maxDeltaRateDegS=max(abs(rad2deg(Y.deltaDot)),[],1);
    Q.allocatorActiveSamples=sum(any(Y.active>0,2));

    Q.maxDTau=max(dTau);
    Q.lateDTau=max(dTau(late));
    Q.maxSurfaceForce=max(dF);

    Q.maxSlideIdentity=max(Y.slidingIdentityResidual);
    Q.maxWindKinematic=max(Y.windKinematicResidual);
    Q.maxQNormError=max(abs(Y.qNorm-1));
    Q.maxQrNormError=max(abs(Y.qrNorm-1));

    % Hard feasibility/invariance criteria.
    dmin=rad2deg(C.P7.deltaMin(:)).';
    dmax=rad2deg(C.P7.deltaMax(:)).';
    rmax=rad2deg(C.P7.deltaRateMax(:)).';

    Q.flagGovernor=(Q.minPhi>=-R.phiTol);
    Q.flagAirspeed=(Q.Vmin>=C.R.Vmin-R.stateTol) && ...
                   (Q.Vmax<=C.R.Vmax+R.stateTol);
    Q.flagGamma=(Q.gammaAbsMaxDeg<=rad2deg(C.R.gammaMax)+1e-8);
    Q.flagAlpha=(Q.alphaMinDeg>=rad2deg(C.R.alphaPhysicalMin)-1e-8) && ...
                (Q.alphaMaxDeg<=rad2deg(C.R.alphaPhysicalMax)+1e-8);
    Q.flagSurfaceBox=all(max(rad2deg(Y.deltaAct),[],1)<=dmax+rad2deg(R.boxTol)) && ...
                     all(min(rad2deg(Y.deltaAct),[],1)>=dmin-rad2deg(R.boxTol));
    Q.flagSurfaceRate=all(Q.maxDeltaRateDegS<=rmax+rad2deg(R.rateTol));
    Q.flagThrustBox=(Q.Tmin>=C.P7.Tmin-R.boxTol) && ...
                    (Q.Tmax<=C.P7.Tmax+R.boxTol);
    Q.flagThrustRate=(Q.maxThrustRate<=C.P7.TRateMax+R.rateTol);
    Q.flagQuaternion=(Q.maxQNormError<=R.qNormTol) && ...
                     (Q.maxQrNormError<=R.qNormTol);
    Q.flagKappaPositive=(Q.kappaMin>0);

    Q.hardFeasible=Q.solverOK && Q.flagGovernor && Q.flagAirspeed && ...
        Q.flagGamma && Q.flagAlpha && Q.flagSurfaceBox && ...
        Q.flagSurfaceRate && Q.flagThrustBox && Q.flagThrustRate && ...
        Q.flagQuaternion && Q.flagKappaPositive;

    % Certified-domain preservation is reported separately because Rev.4.1
    % is an extended robustness study rather than a theorem-premise test.
    Q.flagAttitudeTube=(deg2rad(Q.maxTubeDeg)<=R.totalAttitudeTube+R.stateTol);
    Q.flagBetaTube=(deg2rad(Q.betaAbsMaxDeg)<=R.betaTube+R.stateTol);
    Q.flagKappaCert=(Q.kappaMin>=R.kappaCert-R.stateTol);
    Q.certifiedDomain=Q.hardFeasible && Q.flagAttitudeTube && ...
        Q.flagBetaTube && Q.flagKappaCert;

catch ME
    Q.solverOK=false;
    Q.hardFeasible=false;
    Q.certifiedDomain=false;
    Q.errorId=ME.identifier;
    Q.errorMessage=ME.message;

    if strcmp(ME.identifier,'rg_scalar_governor:backup') || ...
            strcmp(ME.identifier,'rg_raw_command_chart:singular')
        Q.domainExit=true;
        Q.numericalFailure=false;
    else
        Q.domainExit=false;
        Q.numericalFailure=true;
    end
end
end

function Q=local_empty_run()
Q=struct( ...
    'index',0,'solverOK',false,'hardFeasible',false,'certifiedDomain',false, ...
    'domainExit',false,'numericalFailure',false, ...
    'lateErLinf',[nan nan nan],'latePosNormMax',nan,'latePosNormRMS',nan, ...
    'maxPosNorm',nan,'lateWindErrorMax',nan,'minLambda',nan,'minPhi',nan, ...
    'maxTubeDeg',nan,'Vmin',nan,'Vmax',nan,'gammaAbsMaxDeg',nan, ...
    'alphaMinDeg',nan,'alphaMaxDeg',nan,'betaAbsMaxDeg',nan,'kappaMin',nan, ...
    'Tmin',nan,'Tmax',nan,'maxThrustRate',nan, ...
    'maxDeltaActDeg',[nan nan nan],'maxDeltaRateDegS',[nan nan nan], ...
    'allocatorActiveSamples',nan,'maxDTau',nan,'lateDTau',nan, ...
    'maxSurfaceForce',nan,'maxSlideIdentity',nan,'maxWindKinematic',nan, ...
    'maxQNormError',nan,'maxQrNormError',nan, ...
    'flagGovernor',false,'flagAirspeed',false,'flagGamma',false, ...
    'flagAlpha',false,'flagSurfaceBox',false,'flagSurfaceRate',false, ...
    'flagThrustBox',false,'flagThrustRate',false,'flagQuaternion',false, ...
    'flagKappaPositive',false,'flagAttitudeTube',false,'flagBetaTube',false, ...
    'flagKappaCert',false,'errorId','','errorMessage','');
end
