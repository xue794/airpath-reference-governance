function Y=w2c_outputs(t,x,mode,C,method,O)
%W2C_OUTPUTS Common and method-specific diagnostics on the stored output grid.
n=numel(t); nr=18;
Y.er=zeros(n,3);Y.s=zeros(n,3);Y.vc=zeros(n,3);Y.eg=zeros(n,3);Y.ev=zeros(n,3);
Y.tildew=zeros(n,3);Y.what=x(:,17:19);Y.w=zeros(n,3);Y.nw=zeros(n,3);
Y.V=zeros(n,1);Y.gamma=zeros(n,1);Y.chi=zeros(n,1);Y.alpha=zeros(n,1);Y.beta=zeros(n,1);Y.kappa=zeros(n,1);
Y.xf=x(:,20:25);Y.cc=zeros(n,3);Y.u0=zeros(n,3);Y.uGov=zeros(n,3);Y.uB=zeros(n,3);
Y.phi=zeros(n,nr);Y.phi0=zeros(n,nr);Y.correctionNorm=zeros(n,1);Y.objective=zeros(n,1);Y.intervened=false(n,1);
Y.scalarLambda=nan(n,1);Y.scalarLambdaJ=nan(n,nr);Y.scalarLimiter=false(n,nr);
Y.binding=false(n,nr);Y.qpActive=false(n,nr);Y.qpSupport=false(n,nr);Y.qpPath=cell(n,1);
Y.qpPrimal=nan(n,1);Y.qpStationarity=nan(n,1);Y.qpDualViolation=nan(n,1);Y.qpComplementarity=nan(n,1);
Y.qpEnumerationFallback=false(n,1);Y.normalizedPhiViolation=zeros(n,1);
Y.stateMargins=zeros(n,18);
Y.thetaF=zeros(n,1);Y.thetaR=zeros(n,1);Y.thetaC=zeros(n,1);Y.qNorm=zeros(n,1);Y.qrNorm=zeros(n,1);
Y.tauC=zeros(n,3);Y.tauReq=zeros(n,3);Y.dAlloc=zeros(n,3);Y.dServo=zeros(n,3);Y.dTauNom=zeros(n,3);Y.dTauPhysical=zeros(n,3);
Y.deltaCmd=zeros(n,3);Y.deltaAct=x(:,33:35);Y.deltaDot=zeros(n,3);Y.active=zeros(n,3);
Y.Tcmd=zeros(n,1);Y.Tact=x(:,36);Y.TDot=zeros(n,1);Y.fSurface=zeros(n,3);
Y.slidingIdentityResidual=zeros(n,1);Y.windKinematicResidual=zeros(n,1);Y.rd=zeros(n,3);Y.r=x(:,1:3);
Y.constraintNames=w2b_row_names();
Y.stateMarginNames={'hVmin';'hVmax';'hGammaMin';'hGammaMax';'psiVmin';'psiVmax';'psiGammaMin';'psiGammaMax'; ...
 'rateV';'rateGamma';'rateChi';'aGamma';'bank';'Nminus';'Nplus';'Tlower';'Tupper';'backup'};
for k=1:n
    tt=t(k);xx=x(k,:).';S=w2c_control_snapshot(tt,xx,mode,C,method,O);
    U=cl6_plant_input(xx(33:35),xx(36),mode,C);
    xpDot=t2_gna_6dof_quat(tt,xx(1:13),U,C.P);
    Y.er(k,:)=S.Gd.er.';Y.s(k,:)=S.Gd.s.';Y.vc(k,:)=S.Gd.vc.';Y.eg(k,:)=S.eg.';Y.ev(k,:)=S.ev.';
    Y.tildew(k,:)=S.tildew.';Y.w(k,:)=S.Air.wI.';Y.nw(k,:)=S.Mwind.nw.';
    Y.V(k)=S.Air.V;Y.gamma(k)=S.Air.gamma;Y.chi(k)=S.Air.chi;Y.alpha(k)=S.Air.alpha;Y.beta(k)=S.Air.beta;Y.kappa(k)=S.Air.kappa;
    Y.cc(k,:)=S.cc.';Y.u0(k,:)=S.u0.';Y.uGov(k,:)=S.Gov.u.';Y.uB(k,:)=S.Gov.uB.';
    Y.phi(k,:)=S.Gov.phi.';Y.phi0(k,:)=S.Gov.phi0.';Y.correctionNorm(k)=S.Gov.CorrectionNorm;Y.objective(k)=S.Gov.Objective;Y.intervened(k)=S.Gov.Intervened;
    Ck=S.Gov.C; denom=1+abs(Ck.c)+abs(Ck.G)*abs(S.Gov.u);
    Y.normalizedPhiViolation(k)=max([0;-S.Gov.phi./denom]);
    btol=O.BindingRelTol*denom;Y.binding(k,:)=abs(S.Gov.phi).'<=(btol.');
    if strcmpi(method,'Proposed')
        Y.scalarLambda(k)=S.Gov.lambda;Y.scalarLambdaJ(k,:)=S.Gov.lambdaJ.';
        if S.Gov.lambda<1-O.InterventionTol
            lm=min(S.Gov.lambdaJ);lt=O.LimiterTol*(1+abs(lm));Y.scalarLimiter(k,:)=abs(S.Gov.lambdaJ-lm).'<=(lt);
        end
        Y.qpPath{k}='NA';
    else
        D=S.Gov.Solver;Y.qpPath{k}=D.Path;Y.qpEnumerationFallback(k)=D.FallbackUsed;
        if ~isempty(D.ActiveRows),Y.qpActive(k,D.ActiveRows)=true;end
        if ~isempty(D.SupportRows),Y.qpSupport(k,D.SupportRows)=true;end
        Y.qpPrimal(k)=D.PrimalNormalized;Y.qpStationarity(k)=D.StationarityNormalized;
        Y.qpDualViolation(k)=D.DualViolationNormalized;Y.qpComplementarity(k)=D.ComplementarityNormalized;
    end
    nu=xx(23:25); % x_f(4:6) within complete state x(20:25)
    Y.stateMargins(k,:)=[Ck.h(1),Ck.h(2),Ck.h(3),Ck.h(4),Ck.psi1(1),Ck.psi1(2),Ck.psi1(3),Ck.psi1(4), ...
        C.R.nuMax(1)-abs(nu(1)),C.R.nuMax(2)-abs(nu(2)),C.R.nuMax(3)-abs(nu(3)), ...
        Ck.aGamma-C.R.aGammaMin,C.R.muMax-abs(Ck.mu),Ck.hNminus,Ck.hNplus, ...
        Ck.T-(C.R.Tmin+C.R.DeltaT),(C.R.Tmax-C.R.DeltaT)-Ck.T,Ck.hB];
    Y.thetaF(k)=dcm_distance(S.Cc,S.Cr);CIB=S.Air.CIB;Y.thetaR(k)=dcm_distance(S.Cr,CIB);Y.thetaC(k)=dcm_distance(S.Cc,CIB);
    Y.qNorm(k)=norm(xx(7:10));Y.qrNorm(k)=norm(xx(26:29));
    Y.tauC(k,:)=S.tauC.';Y.tauReq(k,:)=S.tauReq.';Y.dAlloc(k,:)=S.dAlloc.';Y.dServo(k,:)=S.dServo.';Y.dTauNom(k,:)=S.dTauNom.';
    Y.deltaCmd(k,:)=S.deltaCmd.';Y.deltaDot(k,:)=S.deltaDot.';Y.active(k,:)=S.Alloc.active.';Y.Tcmd(k)=S.Tcmd;Y.TDot(k)=S.TDot;
    [fA,tauA]=t2_gna_aero(S.Air.vgB,S.Air.omega,xx(33:35),CIB,S.Air.wI,S.Air.h,C.P);
    [f0,~]=t2_gna_aero(S.Air.vgB,S.Air.omega,zeros(3,1),CIB,S.Air.wI,S.Air.h,C.P);
    Y.dTauPhysical(k,:)=(tauA-(-S.Mom.D*S.Air.omega+S.tauC)).';Y.fSurface(k,:)=(fA-f0).';
    sDotActual=xpDot(1:3)-S.Gd.rdDot+C.G.Lambda*S.zeDot;
    sDotTheory=-C.G.Ks*S.Gd.s-C.G.Kc.*S.Gd.sigma+S.tildew+S.ev;
    Y.slidingIdentityResidual(k)=norm(sDotActual-sDotTheory);
    Y.windKinematicResidual(k)=norm(S.Air.vgI-S.Air.vaI-S.Air.wI);Y.rd(k,:)=S.Gd.rd.';
end
end
