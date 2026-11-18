function R8=viii_residual_decomposition_audit(S,C)
%VIII_RESIDUAL_DECOMPOSITION_AUDIT Exact nonoverlapping physical residual split.
%
% IMPORTANT:
% The frozen VI-D/VII nominal maps use rho_ref=C.R.rho, whereas the actual
% nonlinear plant uses ISA density rho(h).  Therefore
%
%   d_tau,physical
%     = B_ref e_delta + d_alloc + d_tau,m,rho,
%
% where d_tau,m,rho is a legitimate model-residual term.  It must NOT be
% forced to zero by the Section-VIII audit.
%
% The same density residual is retained for the surface translational force.

n=numel(S.t);

rotClosureErr=0;
forceClosureErr=0;

maxModelTau=0;
maxModelForce=0;
maxEA=0;
maxDTau=0;
maxFsurfPhysical=0;

for k=1:n
    tt=S.t(k);
    xx=S.x(k,:).';
    Snap=cl6_control_snapshot(tt,xx,S.mode,C);

    deltaAct=xx(33:35);
    eDelta=deltaAct-Snap.deltaCmd;
    eT=xx(36)-Snap.Tcmd;

    % Internal actuator/allocation contribution in the frozen nominal map.
    dTauInternal=Snap.Mom.B*eDelta+Snap.dAlloc;

    % Exact density-scheduling/model residual:
    % physical GNA uses rho(h), nominal allocator uses rho_ref.
    rhoAct=t2_isa_density(Snap.Air.h);
    rhoRef=C.R.rho;
    rhoRatio=rhoAct/rhoRef;

    tauNomFull=Snap.Mom.tau0 + Snap.Mom.K*Snap.Air.omega + ...
               Snap.Mom.B*deltaAct;
    dTauModel=(rhoRatio-1)*tauNomFull;

    dTauPhysical=S.Y.dTauPhysical(k,:).';
    rotClosureErr=max(rotClosureErr, ...
        norm(dTauPhysical-(dTauInternal+dTauModel)));
    maxModelTau=max(maxModelTau,norm(dTauModel));

    % Translational control-surface force:
    % command-level + actuator-lag + density-model residual.
    fCmdNom=vii_surface_force_residual( ...
        Snap.Air.V,Snap.Air.alpha,Snap.deltaCmd,C.R);
    fActNom=vii_surface_force_residual( ...
        Snap.Air.V,Snap.Air.alpha,deltaAct,C.R);
    fLagNom=fActNom-fCmdNom;
    dFModel=(rhoRatio-1)*fActNom;

    fPhysical=S.Y.fSurface(k,:).';
    forceClosureErr=max(forceClosureErr, ...
        norm(fPhysical-(fCmdNom+fLagNom+dFModel)));
    maxModelForce=max(maxModelForce,norm(dFModel));

    maxEA=max(maxEA,norm([eDelta;eT]));
    maxDTau=max(maxDTau,norm(dTauPhysical));
    maxFsurfPhysical=max(maxFsurfPhysical,norm(fPhysical));
end

% Backward-compatible field names rotErr/forceErr now mean the FULL
% nonoverlapping closure residual, including the explicit model term.
R8.rotErr=rotClosureErr;
R8.forceErr=forceClosureErr;
R8.maxModelTau=maxModelTau;
R8.maxModelForce=maxModelForce;
R8.maxEA=maxEA;
R8.maxDTau=maxDTau;
R8.maxFsurf=maxFsurfPhysical;
end
