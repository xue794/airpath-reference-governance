function E=rev4_extended_metrics(CL,Cert)
%REV4_EXTENDED_METRICS Consolidated paper-facing complete-6DOF metrics.
%
% Only certified_small is used as the local theorem-instance case.
% The other scenarios are extended nonlinear robustness/performance tests.

C=CL.C;

% MATLAB R2018a requires the destination struct array to have the same
% field set before indexed assignment.  Initialize the array with the first
% fully populated case rather than with struct([]).
Case=local_case_metrics(CL.hold,C,'hold', ...
    'baseline-regression','still-air hold / implementation sanity');
Case(2)=local_case_metrics(CL.constantCold,C,'constant_cold', ...
    'extended-robustness','constant-wind observer cold-start transient');
Case(3)=local_case_metrics(CL.varying,C,'varying', ...
    'extended-robustness','time-varying wind, full nonlinear loop');
Case(4)=local_case_metrics(CL.noisy,C,'noisy', ...
    'extended-robustness','time-varying wind + bounded measurement error');
Case(5)=local_case_metrics(CL.wrap,C,'wrap', ...
    'extended-robustness','S^1 course branch crossing');
Case(6)=local_case_metrics(Cert.Sim,Cert.C,'certified_small', ...
    'theorem-instance','strict-interior local macro-cascade numerical instance');

E=struct();
E.cases=Case;

% Dedicated varying-wind physical residual summary.
V=CL.varying;
iv=V.t>=V.Case.tSteady;
dTau=sqrt(sum(V.Y.dTauPhysical.^2,2));
dAlloc=sqrt(sum(V.Y.dAlloc.^2,2));
fSurf=sqrt(sum(V.Y.fSurface.^2,2));

E.varyingPhysical=struct( ...
    'maxDTau',max(dTau), ...
    'lateDTau',max(dTau(iv)), ...
    'maxDAlloc',max(dAlloc), ...
    'maxSurfaceForce',max(fSurf), ...
    'maxThrustLag',max(abs(V.Y.Tact-V.Y.Tcmd)), ...
    'maxThrustRate',max(abs(V.Y.TDot)), ...
    'allocatorActiveSamples',sum(any(V.Y.active>0,2)));

fprintf('\nREV.4 CONSOLIDATED COMPLETE-6DOF METRICS\n');
fprintf('========================================\n');
for k=1:numel(Case)
    Q=Case(k);
    fprintf('%-16s %-19s late L_inf e_r=[%.4g %.4g %.4g] m, ', ...
        Q.mode,Q.role,Q.latePosLinf);
    fprintf('tube=%.5f deg, minPhi=%.3e, minLambda=%.4f\n', ...
        Q.maxTubeDeg,Q.minPhi,Q.minLambda);
end

fprintf('varying physical max/late ||d_tau|| = [%.6f %.6f] N m\n', ...
    E.varyingPhysical.maxDTau,E.varyingPhysical.lateDTau);
fprintf('varying max omitted surface-force residual = %.6f N\n', ...
    E.varyingPhysical.maxSurfaceForce);
end

function Q=local_case_metrics(S,C,mode,role,description)
late=S.t>=S.Case.tSteady;
if ~any(late), late=true(size(S.t)); end

er=S.Y.er;
posNorm=sqrt(sum(er.^2,2));
windErr=sqrt(sum(S.Y.tildew.^2,2));
dTau=sqrt(sum(S.Y.dTauPhysical.^2,2));
dAlloc=sqrt(sum(S.Y.dAlloc.^2,2));
fSurf=sqrt(sum(S.Y.fSurface.^2,2));

Q=struct();
Q.mode=mode;
Q.role=role;
Q.description=description;

Q.latePosLinf=max(abs(er(late,:)),[],1);
Q.latePosRMS=sqrt(mean(er(late,:).^2,1));
Q.latePosNormMax=max(posNorm(late));
Q.lateWindErrorMax=max(windErr(late));

Q.minLambda=min(S.Y.lambda);
Q.minPhi=min(S.Y.minPhi);

Q.maxThetaFDeg=rad2deg(max(S.Y.thetaF));
Q.maxThetaRDeg=rad2deg(max(S.Y.thetaR));
Q.maxTubeDeg=rad2deg(max(S.Y.thetaF+S.Y.thetaR));

Q.Vmin=min(S.Y.V);
Q.Vmax=max(S.Y.V);
Q.gammaAbsMaxDeg=rad2deg(max(abs(S.Y.gamma)));
Q.alphaMinDeg=rad2deg(min(S.Y.alpha));
Q.alphaMaxDeg=rad2deg(max(S.Y.alpha));
Q.betaAbsMaxDeg=rad2deg(max(abs(S.Y.beta)));
Q.kappaMin=min(S.Y.kappa);

Q.Tmin=min(S.Y.Tact);
Q.Tmax=max(S.Y.Tact);
Q.maxThrustLag=max(abs(S.Y.Tact-S.Y.Tcmd));
Q.maxThrustRate=max(abs(S.Y.TDot));

Q.maxDeltaActDeg=max(abs(rad2deg(S.Y.deltaAct)),[],1);
Q.maxDeltaRateDegS=max(abs(rad2deg(S.Y.deltaDot)),[],1);
Q.allocatorActiveSamples=sum(any(S.Y.active>0,2));

Q.maxDTau=max(dTau);
Q.lateDTau=max(dTau(late));
Q.maxDAlloc=max(dAlloc);
Q.maxSurfaceForce=max(fSurf);

Q.maxSlidingIdentityResidual=max(S.Y.slidingIdentityResidual);
Q.maxWindKinematicResidual=max(S.Y.windKinematicResidual);
Q.maxPlantQuatNormError=max(abs(S.Y.qNorm-1));
Q.maxFilterQuatNormError=max(abs(S.Y.qrNorm-1));

% Late-window SMC position-bound proxy.  This is not presented as a global
% continuous-domain theorem bound.
rs=max(abs(S.Y.s(late,:)),[],1);
rg=max(abs(S.Y.eg(late,:)),[],1);
lam=diag(C.G.Lambda).';
Q.positionBoundProxy=2*rs+rg./lam;
Q.positionProxyRatio=Q.latePosLinf./Q.positionBoundProxy;
end
