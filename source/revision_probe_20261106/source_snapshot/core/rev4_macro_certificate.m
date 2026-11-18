function Cert=rev4_macro_certificate()
%REV4_MACRO_CERTIFICATE Final sampled T-2 instantiation of the analytic
% local macro-cascade stability certificate.
%
% This routine does NOT define the theory.  It numerically instantiates the
% repaired theory on the frozen T-2 GTM--GNA local validation case.
%
% x_C = col(z,s,xi_p,x_G,x_F,x_R,x_A),  dim=31
% x_w = tilde w,                         dim=3
%
% The estimator-error dynamics are autonomous analytically, so the final
% two-block comparison graph is a cascade:
%       x_w ---> x_C
%       x_C -/-> x_w
% and rho(D_M^{-1}G_M)=0 exactly.

R=rev4_params();
C=cl6_params();
C.CertScale=R.certScale;
P8=viii_params();
B=viii_blocks();

idxC=[];
for k=1:numel(R.macroPrimitiveIDs)
    idxC=[idxC B.idx{R.macroPrimitiveIDs(k)}]; %#ok<AGROW>
end
idxW=B.idx{R.windPrimitiveID};

SC=diag(P8.scaleVector(idxC));
SW=diag(P8.scaleVector(idxW));

X0=zeros(B.n,1);

fprintf('\nREV.4 LOCAL MACRO-CASCADE CERTIFICATE\n');
fprintf('=====================================\n');
fprintf('certified-small scale = %.8f\n',R.certScale);
fprintf('macro/wind dimensions = %d / %d\n',numel(idxC),numel(idxW));

% Run the local full-6DOF case once.
S=simulate_complete_6dof_case(R.certMode,C);

Fam=struct('A',{},'kind',{},'t',{},'X',{}, ...
           'maxReal',{},'q',{},'lambda',{},'minPhi',{}, ...
           'normX',{},'allocatorActive',{},'tube',{});

% Nominal strict-interior hold point.
A=local_block(R.certHoldTime,X0,'hold',idxC,idxC,R.relStep,C,P8);
Fam(end+1)=local_family_entry(SC\A*SC,'hold',R.certHoldTime,X0, ...
    nan,nan,0,false,0); %#ok<AGROW>

% Zero-error reference-lifting family.
refMinLambda=inf;
refMinPhi=inf;

for kt=1:numel(R.certReferenceTimes)
    t=R.certReferenceTimes(kt);

    [xf,~]=viii_lift_reduced(t,X0,R.certMode,C,P8);
    Xback=viii_reduced_from_full(t,xf,R.certMode,C,P8);

    rt=norm((Xback-X0)./P8.scaleVector,inf);
    if rt>R.roundTripTol
        error('Rev4:Lift','Reference-lifting round-trip error %.3e at t=%g.',rt,t);
    end

    Snap=cl6_control_snapshot(t,xf,R.certMode,C);
    lam=Snap.Gov.lambda;
    phi=min(Snap.Gov.phi);

    if phi<-R.phiTol
        error('Rev4:ReferenceFeasibility', ...
            'Certified-small reference lifting lost feasibility at t=%g, min phi=%.6e.',t,phi);
    end

    refMinLambda=min(refMinLambda,lam);
    refMinPhi=min(refMinPhi,phi);

    A=local_block(t,X0,R.certMode,idxC,idxC,R.relStep,C,P8);
    Fam(end+1)=local_family_entry(SC\A*SC,'reference',t,X0, ...
        lam,phi,0,false,0); %#ok<AGROW>
end

% Realized local full-loop family.
realMinLambda=inf;
realMinPhi=inf;
realMaxNormX=0;
realAllocatorActive=0;
realMaxTube=0;

for kt=1:numel(R.certRealizedTimes)
    [~,ii]=min(abs(S.t-R.certRealizedTimes(kt)));
    t=S.t(ii);
    X=viii_reduced_from_full(t,S.x(ii,:).',R.certMode,C,P8);

    [xf,~]=viii_lift_reduced(t,X,R.certMode,C,P8);
    Snap=cl6_control_snapshot(t,xf,R.certMode,C);

    lam=Snap.Gov.lambda;
    phi=min(Snap.Gov.phi);
    nX=norm(X./P8.scaleVector,2);
    act=any(Snap.Alloc.active);

    thetaF=dcm_distance(Snap.Cc,Snap.Cr);
    thetaR=dcm_distance(Snap.Cr,Snap.Air.CIB);
    tube=thetaF+thetaR;

    realMinLambda=min(realMinLambda,lam);
    realMinPhi=min(realMinPhi,phi);
    realMaxNormX=max(realMaxNormX,nX);
    realAllocatorActive=realAllocatorActive+act;
    realMaxTube=max(realMaxTube,tube);

    A=local_block(t,X,R.certMode,idxC,idxC,R.relStep,C,P8);
    Fam(end+1)=local_family_entry(SC\A*SC,'realized',t,X, ...
        lam,phi,nX,act,tube); %#ok<AGROW>
end

N=numel(Fam);
if N~=R.expectedFamilyCount
    error('Rev4:FamilyCount','Expected %d macro matrices, obtained %d.', ...
        R.expectedFamilyCount,N);
end

maxReal=-inf;
for k=1:N
    Fam(k).maxReal=max(real(eig(Fam(k).A)));
    maxReal=max(maxReal,Fam(k).maxReal);
end

fprintf('MC-P1 sampled macro family:\n');
fprintf('   total matrices = %d\n',N);
fprintf('   hold/reference/realized = %d / %d / %d\n', ...
    sum(strcmp({Fam.kind},'hold')), ...
    sum(strcmp({Fam.kind},'reference')), ...
    sum(strcmp({Fam.kind},'realized')));
fprintf('   worst max Re eig(A_C,k) = %.9e\n',maxReal);
fprintf('   reference min lambda / min phi = [%.9f %.6e]\n',refMinLambda,refMinPhi);
fprintf('   realized min lambda / min phi = [%.9f %.6e]\n',realMinLambda,realMinPhi);
fprintf('   realized max normalized ||X_red|| = %.6f\n',realMaxNormX);
fprintf('   realized allocator-active sampled centers = %d / %d\n', ...
    realAllocatorActive,numel(R.certRealizedTimes));
fprintf('   realized max total attitude tube = %.6f deg\n',rad2deg(realMaxTube));

assert(maxReal<-R.hurwitzTol,'Rev4 macro family contains a non-Hurwitz sampled matrix.');
assert(refMinPhi>-R.phiTol && realMinPhi>-R.phiTol,'Rev4 strict-interior feasibility check failed.');
assert(realAllocatorActive==0,'Rev4 local theorem-instance unexpectedly activated the surface allocator.');

% ---------------------------------------------------------------------
% Common quadratic macro storage LMI.
% ---------------------------------------------------------------------
n=numel(idxC);

setlmis([]);
Pvar=lmivar(1,[n 1]);

% I - P < 0.
ilmi=1;
lmiterm([ilmi 1 1 Pvar],-1,1);
lmiterm([ilmi 1 1 0],eye(n));

% A_k' P + P A_k + I < 0.
for k=1:N
    ilmi=ilmi+1;
    Ak=Fam(k).A;
    lmiterm([ilmi 1 1 Pvar],Ak.',1,'s');
    lmiterm([ilmi 1 1 0],eye(n));
end

lmis=getlmis;
[tmin,xfeas]=feasp(lmis);

fprintf('MC-P2 LMI common-storage search:\n');
fprintf('   feasp tmin = %.9e\n',tmin);

if isempty(xfeas) || ~isfinite(tmin) || tmin>=0
    error('Rev4:CommonP','Common quadratic macro storage was not found.');
end

Praw=dec2mat(lmis,xfeas,Pvar);
Praw=0.5*(Praw+Praw.');

% Normalize only for reporting/export.  Positivity and dissipation signs are
% invariant under positive scalar scaling.
P=Praw*(n/trace(Praw));

qmin=inf;
kworst=1;
for k=1:N
    Q=-(Fam(k).A.'*P+P*Fam(k).A);
    Q=0.5*(Q+Q.');
    Fam(k).q=min(eig(Q));
    if Fam(k).q<qmin
        qmin=Fam(k).q;
        kworst=k;
    end
end

minP=min(eig(P));
maxP=max(eig(P));
condP=cond(P);
decayProxy=qmin/maxP;

fprintf('MC-P3 direct normalized common-P verification:\n');
fprintf('   trace(P_C) = %.6f (normalized to macro dimension)\n',trace(P));
fprintf('   min eig(P_C) / cond(P_C) = [%.6e %.6e]\n',minP,condP);
fprintf('   worst q_k = %.9e\n',qmin);
fprintf('   worst point = %s, t=%.3f s\n',Fam(kworst).kind,Fam(kworst).t);
fprintf('   scale-invariant q_min/lambda_max(P_C) = %.9e 1/s\n',decayProxy);

assert(minP>0 && qmin>R.commonQTol,'Rev4 direct common-P verification failed.');

% Worst-point finite-difference convergence.
qStep=nan(size(R.stepCheck));
for ks=1:numel(R.stepCheck)
    Fw=Fam(kworst);

    if strcmp(Fw.kind,'hold')
        mode='hold';
    else
        mode=R.certMode;
    end

    Aw=local_block(Fw.t,Fw.X,mode,idxC,idxC,R.stepCheck(ks),C,P8);
    Aws=SC\Aw*SC;
    Qw=-(Aws.'*P+P*Aws);
    Qw=0.5*(Qw+Qw.');
    qStep(ks)=min(eig(Qw));
end

fprintf('   worst-point step check q = [');
fprintf(' %.9e',qStep);
fprintf(' ]\n');

assert(all(qStep>R.commonQTol),'Rev4 finite-difference convergence check lost common-P positivity.');

% ---------------------------------------------------------------------
% Causal wind-estimator interconnection synchronization.
% ---------------------------------------------------------------------
maxMacroToWind=0;
maxWindToMacro=0;
maxPCWindCoupling=0;

for k=1:N
    Fk=Fam(k);

    if strcmp(Fk.kind,'hold')
        mode='hold';
    else
        mode=R.certMode;
    end

    AWC=local_block(Fk.t,Fk.X,mode,idxW,idxC,R.relStep,C,P8);
    ACW=local_block(Fk.t,Fk.X,mode,idxC,idxW,R.relStep,C,P8);

    AWCs=SW\AWC*SC;
    ACWs=SC\ACW*SW;

    maxMacroToWind=max(maxMacroToWind,norm(AWCs,2));
    maxWindToMacro=max(maxWindToMacro,norm(ACWs,2));
    maxPCWindCoupling=max(maxPCWindCoupling,2*norm(P*ACWs,2));
end

rhoMacro=0; % Exact analytic result: autonomous estimator error dynamics.

fprintf('MC-P4 macro-cascade synchronization:\n');
fprintf('   max numerical ||d f_w / d x_C||_2 = %.9e\n',maxMacroToWind);
fprintf('   max numerical ||d f_C / d x_w||_2 = %.9e\n',maxWindToMacro);
fprintf('   max 2||P_C A_Cw||_2 = %.9e\n',maxPCWindCoupling);
fprintf('   analytic macro -> wind state gain = 0 exactly\n');
fprintf('   analytic rho(D_M^{-1}G_M) = %.1f\n',rhoMacro);

assert(maxMacroToWind<R.macroToWindTol, ...
    'Rev4 reduced-vector-field implementation violates autonomous estimator-error structure.');

Cert=struct();
Cert.status='REV4_MACRO_CASCADE_CERTIFICATE_PASSED';
Cert.R=R;
Cert.C=C;
Cert.P8=P8;
Cert.Sim=S;
Cert.Family=Fam;
Cert.P=P;
Cert.Praw=Praw;
Cert.tmin=tmin;
Cert.minP=minP;
Cert.condP=condP;
Cert.qmin=qmin;
Cert.decayProxy=decayProxy;
Cert.worstIndex=kworst;
Cert.qStep=qStep;
Cert.maxReal=maxReal;
Cert.refMinLambda=refMinLambda;
Cert.refMinPhi=refMinPhi;
Cert.realMinLambda=realMinLambda;
Cert.realMinPhi=realMinPhi;
Cert.realMaxNormX=realMaxNormX;
Cert.realAllocatorActive=realAllocatorActive;
Cert.realMaxTube=realMaxTube;
Cert.maxMacroToWind=maxMacroToWind;
Cert.maxWindToMacro=maxWindToMacro;
Cert.maxPCWindCoupling=maxPCWindCoupling;
Cert.rhoMacro=rhoMacro;

fprintf('\nREV.4 MACRO-CASCADE CERTIFICATE PASSED.\n');
end

function F=local_family_entry(A,kind,t,X,lambda,minPhi,normX,allocatorActive,tube)
F=struct('A',A,'kind',kind,'t',t,'X',X, ...
         'maxReal',nan,'q',nan,'lambda',lambda,'minPhi',minPhi, ...
         'normX',normX,'allocatorActive',allocatorActive,'tube',tube);
end

function A=local_block(t,X,mode,ii,jj,relStep,C,P8)
A=zeros(numel(ii),numel(jj));

for kc=1:numel(jj)
    k=jj(kc);
    h=relStep*P8.scaleVector(k);
    if h<=0, h=relStep; end

    Xp=X;
    Xm=X;
    Xp(k)=Xp(k)+h;
    Xm(k)=Xm(k)-h;

    fp=viii_reduced_vector_field(t,Xp,mode,C,P8);
    fm=viii_reduced_vector_field(t,Xm,mode,C,P8);

    A(:,kc)=(fp(ii)-fm(ii))/(2*h);
end
end
