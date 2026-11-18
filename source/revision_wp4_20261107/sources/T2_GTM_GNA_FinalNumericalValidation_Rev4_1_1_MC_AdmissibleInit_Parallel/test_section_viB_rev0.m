function Results=test_section_viB_rev0()
%TEST_SECTION_VIB_REV0 Section VI-B geometry tests VIB-P1--VIB-P8.

VIA=test_section_viA_rev0();
R=VIA.R;
A=VIA.A;
B=viB_params();

fprintf('\nSection VI-B Command Attitude Geometry + AoA/Sideslip Transfer Rev.0\n');
fprintf('====================================================================\n');

%% VIB-P1 command attitude geometry
rng(31,'twister');
maxSO3=0; maxAxis=0; maxBc=0;
for k=1:1000
    gamma=deg2rad(-8+16*rand);
    chi=-pi+2*pi*rand;
    mu=deg2rad(-30+60*rand);
    alpha=R.alphaLo+(R.alphaHi-R.alphaLo)*rand;

    P=viA_airpath_basis(gamma,chi);
    Cc=viB_command_attitude(gamma,chi,mu,alpha);
    ell=cos(mu)*P.nGamma+sin(mu)*P.nChi;

    target=cos(alpha)*P.t+sin(alpha)*ell;
    bcTarget=[cos(alpha);0;sin(alpha)];

    maxSO3=max(maxSO3,max(norm(Cc.'*Cc-eye(3),'fro'),abs(det(Cc)-1)));
    maxAxis=max(maxAxis,norm(Cc(:,1)-target));
    maxBc=max(maxBc,norm(Cc.'*P.t-bcTarget));
end
fprintf('VIB-P1 command-attitude identity:\n');
fprintf('   max SO(3) residual          = %.3e\n',maxSO3);
fprintf('   max Cc*e1 axis residual     = %.3e\n',maxAxis);
fprintf('   max b_c identity residual   = %.3e\n',maxBc);
assert(maxSO3<1e-12 && maxAxis<1e-12 && maxBc<1e-12,'VIB-P1 failed');

%% VIB-P2 relative-rotation composition and airdata reconstruction
rng(32,'twister');
maxComp=0; maxB=0; maxInertial=0; maxZeroA=0; maxZeroB=0;
for k=1:1000
    gamma=deg2rad(-5+10*rand); chi=-pi+2*pi*rand;
    mu=deg2rad(-20+40*rand); alphaC=R.alphaLo+(R.alphaHi-R.alphaLo)*rand;

    ax=randn(3,1); ax=ax/norm(ax);
    ay=randn(3,1); ay=ay/norm(ay);
    Qf=rotvec_to_dcm(ax*(B.thetaFMax*rand));
    QR=rotvec_to_dcm(ay*(B.thetaRMax*rand));

    G=viB_geometry(gamma,chi,mu,alphaC,Qf,QR);
    maxComp=max(maxComp,norm(G.Cc.'*G.C-G.Qc,'fro'));

    bRec=[cos(G.alpha)*cos(G.beta);
          sin(G.beta);
          sin(G.alpha)*cos(G.beta)];
    maxB=max(maxB,norm(bRec-G.b));
    maxInertial=max(maxInertial,norm(G.C*bRec-G.B.t));

    G0=viB_geometry(gamma,chi,mu,alphaC,eye(3),eye(3));
    maxZeroA=max(maxZeroA,abs(rg_angle_error(G0.alpha,alphaC)));
    maxZeroB=max(maxZeroB,abs(G0.beta));
end
fprintf('VIB-P2 relative rotations / airdata reconstruction:\n');
fprintf('   max Qc composition residual = %.3e\n',maxComp);
fprintf('   max body-direction residual  = %.3e\n',maxB);
fprintf('   max inertial-direction residual = %.3e\n',maxInertial);
fprintf('   zero-error max [|alpha-alpha_c| |beta|] = [%.3e %.3e] rad\n',maxZeroA,maxZeroB);
assert(max([maxComp maxB maxInertial maxZeroA maxZeroB])<1e-12,'VIB-P2 failed');

%% VIB-P3 local attitude-tube beta/AoA bounds
rng(33,'twister');
maxTri=0; maxBetaRatio=0; maxAlphaRatio=0;
maxBeta=0; maxAlphaErr=0; minKappa=inf;
for k=1:B.Nmc
    gamma=deg2rad(-8+16*rand); chi=-pi+2*pi*rand;
    mu=deg2rad(-30+60*rand); alphaC=R.alphaLo+(R.alphaHi-R.alphaLo)*rand;

    ax=randn(3,1); ax=ax/norm(ax);
    ay=randn(3,1); ay=ay/norm(ay);
    thF=B.thetaFMax*rand; thR=B.thetaRMax*rand;
    G=viB_geometry(gamma,chi,mu,alphaC, ...
        rotvec_to_dcm(ax*thF),rotvec_to_dcm(ay*thR));

    tsum=G.thetaF+G.thetaR;
    maxTri=max(maxTri,G.thetaC-tsum);
    maxBeta=max(maxBeta,abs(G.beta));
    da=abs(rg_angle_error(G.alpha,alphaC));
    maxAlphaErr=max(maxAlphaErr,da);
    minKappa=min(minKappa,G.kappa);

    if tsum>1e-12
        maxBetaRatio=max(maxBetaRatio,abs(G.beta)/tsum);
        maxAlphaRatio=max(maxAlphaRatio,da/(tsum/cos(B.thetaBar)));
    end
end
fprintf('VIB-P3 local attitude-tube transfer bounds:\n');
fprintf('   max(theta_c-theta_f-theta_R) = %.3e rad\n',maxTri);
fprintf('   max |beta| = %.6f deg, bound %.6f deg, ratio %.4f\n', ...
    rad2deg(maxBeta),rad2deg(B.betaBar),maxBetaRatio);
fprintf('   max |alpha-alpha_c| = %.6f deg, bound %.6f deg, ratio %.4f\n', ...
    rad2deg(maxAlphaErr),rad2deg(B.deltaAlphaMax),maxAlphaRatio);
fprintf('   sampled min kappa_T = %.9f\n',minKappa);
assert(maxTri<1e-10 && maxBeta<=B.betaBar+1e-10 && ...
       maxAlphaErr<=B.deltaAlphaMax+1e-10 && ...
       maxBetaRatio<=1+1e-8 && maxAlphaRatio<=1+1e-8,'VIB-P3 failed');

%% VIB-P4 kappa positivity and geometric thrust transfer
rng(34,'twister');
minK=inf; maxTdiff=0; maxBoundRatio=0;
for k=1:5000
    gamma=deg2rad(-8+16*rand); chi=-pi+2*pi*rand;
    mu=deg2rad(-30+60*rand); alphaC=R.alphaLo+(R.alphaHi-R.alphaLo)*rand;

    ax=randn(3,1); ax=ax/norm(ax);
    ay=randn(3,1); ay=ay/norm(ay);
    thF=B.thetaFMax*rand; thR=B.thetaRMax*rand;
    G=viB_geometry(gamma,chi,mu,alphaC, ...
        rotvec_to_dcm(ax*thF),rotvec_to_dcm(ay*thR));

    Fp=10+20*rand;
    T0=Fp/cos(alphaC);
    T=Fp/G.kappa;
    dt=abs(T-T0);
    tsum=G.thetaF+G.thetaR;

    % Uniform local bound:
    % |kappa-kappa0| <= |dalpha| + beta^2/2
    % <= thetaSum/cos(thetaBar) + thetaSum^2/2.
    db=Fp*(tsum/cos(B.thetaBar)+0.5*tsum^2)/(B.kappaMin*cos(R.alphaHi));

    minK=min(minK,G.kappa);
    maxTdiff=max(maxTdiff,dt);
    if db>1e-14
        maxBoundRatio=max(maxBoundRatio,dt/db);
    end
end
fprintf('VIB-P4 body-axis thrust-port geometry:\n');
fprintf('   certified kappa_min = %.9f, sampled min = %.9f\n',B.kappaMin,minK);
fprintf('   max geometric |T-T_ideal| = %.6f N\n',maxTdiff);
fprintf('   max realized/conservative mismatch-bound ratio = %.4f\n',maxBoundRatio);
assert(minK>=B.kappaMin-1e-10 && maxBoundRatio<=1+1e-8,'VIB-P4 failed');

%% VIB-P5 sideslip derivative identity and beta=0 non-invariance
rng(35,'twister');
maxBetaDotErr=0;
for k=1:B.Nder
    gamma=deg2rad(-6+12*rand); chi=-pi+2*pi*rand;
    mu=deg2rad(-25+50*rand); alphaC=deg2rad(2+5*rand);
    ax=randn(3,1); ax=ax/norm(ax);
    Q=rotvec_to_dcm(ax*(B.thetaTube*(2*rand-1)));
    G=viB_geometry(gamma,chi,mu,alphaC,Q,eye(3));

    omega=(-0.4+0.8*rand(3,1));
    gammaDot=deg2rad(-5+10*rand);
    chiDot=deg2rad(-8+16*rand);
    tdot=viB_tangent_rate(gamma,chi,gammaDot,chiDot);

    bd=viB_beta_rate(G.C,G.B.t,omega,tdot);

    dt=B.fdDt;
    C1=G.C*rotvec_to_dcm(omega*dt);
    P1=viA_airpath_basis(gamma+gammaDot*dt,chi+chiDot*dt);
    A1=viB_airdata_from_geometry(C1,P1.t);
    bfd=(A1.beta-G.beta)/dt;
    maxBetaDotErr=max(maxBetaDotErr,abs(bfd-bd));
end

gamma=0; chi=0; mu=0; alphaC=deg2rad(5);
G0=viB_geometry(gamma,chi,mu,alphaC,eye(3),eye(3));
omega0=[0;0;0.2];
bd0=viB_beta_rate(G0.C,G0.B.t,omega0,zeros(3,1));

fprintf('VIB-P5 sideslip derivative / non-invariance:\n');
fprintf('   max finite-difference beta_dot error = %.3e rad/s\n',maxBetaDotErr);
fprintf('   beta=0 example beta_dot = %.6f rad/s (must be nonzero)\n',bd0);
assert(maxBetaDotErr<2e-6 && abs(bd0)>0.1,'VIB-P5 failed');

%% VIB-P6 dynamic synthetic attitude-error feasibility transfer
S6=simulate_viA_case('aggressive',A.Taggressive,R,A,wind_estimator_params());
n=numel(S6.t);
alpha=zeros(n,1); beta=zeros(n,1); kappa=zeros(n,1); Tphys=zeros(n,1);
thetaF=zeros(n,1); thetaR=zeros(n,1); thetaC=zeros(n,1);
for k=1:n
    y=S6.y(k,:).';
    Syn=viB_synthetic_relative_rotations(S6.t(k),B);
    G=viB_geometry(y(2),y(3),S6.mu(k),S6.alpha(k),Syn.Qf,Syn.QR);
    alpha(k)=G.alpha; beta(k)=G.beta; kappa(k)=G.kappa;
    Tphys(k)=S6.Fparallel(k)/G.kappa;
    thetaF(k)=G.thetaF; thetaR(k)=G.thetaR; thetaC(k)=G.thetaC;
end

amax=max(alpha); amin=min(alpha); bmax=max(abs(beta));
Tmax=max(Tphys); Tmin=min(Tphys); kmin=min(kappa);
dA=max(abs(arrayfun(@(k)rg_angle_error(alpha(k),S6.alpha(k)),(1:n).')));
dT=max(abs(Tphys-S6.Tstar));
maxTube=max(thetaF+thetaR);

fprintf('VIB-P6 aggressive-reference geometric feasibility transfer:\n');
fprintf('   max theta_f+theta_R = %.6f deg (tube %.6f deg)\n', ...
    rad2deg(maxTube),rad2deg(B.thetaTube));
fprintf('   actual alpha range  = [%.6f %.6f] deg\n',rad2deg(amin),rad2deg(amax));
fprintf('   max |beta|          = %.6f deg\n',rad2deg(bmax));
fprintf('   min kappa_T         = %.9f\n',kmin);
fprintf('   physical T* range   = [%.6f %.6f] N\n',Tmin,Tmax);
fprintf('   max geometry-induced [|dalpha| |dT|] = [%.6f deg %.6f N]\n', ...
    rad2deg(dA),dT);

assert(maxTube<=B.thetaTube+1e-10 && ...
       amin>R.alphaPhysicalMin && amax<R.alphaPhysicalMax && ...
       bmax<=B.betaBar+1e-10 && kmin>=B.kappaMin-1e-10 && ...
       Tmin>R.Tmin+R.DeltaT && Tmax<R.Tmax-R.DeltaT,'VIB-P6 failed');

%% VIB-P7 S1 branch invariance of command attitude
S7=simulate_viA_case('wrap',A.Twrap,R,A,wind_estimator_params());
n=numel(S7.t);
maxWrapC=0; maxStepC=0; rawJump=max(abs(diff(S7.cc(:,3))));
Cprev=[];
for k=1:n
    y=S7.y(k,:).';
    C1=viB_command_attitude(y(2),y(3),S7.mu(k),S7.alpha(k));
    chiW=atan2(sin(y(3)),cos(y(3)));
    C2=viB_command_attitude(y(2),chiW,S7.mu(k),S7.alpha(k));
    maxWrapC=max(maxWrapC,dcm_distance(C1,C2));
    if ~isempty(Cprev)
        maxStepC=max(maxStepC,dcm_distance(Cprev,C1));
    end
    Cprev=C1;
end
fprintf('VIB-P7 S1 command-attitude branch invariance:\n');
fprintf('   raw atan2 command-course jump = %.4f deg\n',rad2deg(rawJump));
fprintf('   max Cc(chi)-Cc(wrap(chi)) geodesic = %.3e rad\n',maxWrapC);
fprintf('   max consecutive command-attitude step = %.6f deg\n',rad2deg(maxStepC));
assert(maxWrapC<5e-8 && maxStepC<deg2rad(0.2),'VIB-P7 failed');

%% VIB-P8 deterministic certified-tube boundary stress audit
axes=[eye(3) -eye(3)];
minK8=inf; maxA8=-inf; minA8=inf; maxB8=0; maxT8=0;
for ia=1:size(axes,2)
    for ib=1:size(axes,2)
        for sa=[-1 1]
            for sb=[-1 1]
                Qf=rotvec_to_dcm(axes(:,ia)*(sa*B.thetaFMax));
                QR=rotvec_to_dcm(axes(:,ib)*(sb*B.thetaRMax));
                for alphaC=[R.alphaLo R.alphaHi]
                    G=viB_geometry(deg2rad(4),deg2rad(35),deg2rad(20),alphaC,Qf,QR);
                    minK8=min(minK8,G.kappa);
                    maxA8=max(maxA8,G.alpha); minA8=min(minA8,G.alpha);
                    maxB8=max(maxB8,abs(G.beta));

                    % Stress the largest Section-V retained nominal thrust margin.
                    Tideal=R.Tmax-R.DeltaT;
                    Fp=Tideal*cos(alphaC);
                    Tactual=Fp/G.kappa;
                    maxT8=max(maxT8,Tactual);
                end
            end
        end
    end
end
fprintf('VIB-P8 certified-tube boundary stress audit:\n');
fprintf('   alpha actual extrema = [%.6f %.6f] deg\n',rad2deg(minA8),rad2deg(maxA8));
fprintf('   max |beta|           = %.6f deg\n',rad2deg(maxB8));
fprintf('   min kappa_T          = %.9f (certificate %.9f)\n',minK8,B.kappaMin);
fprintf('   max actual thrust from T_ff=57 N stress = %.6f N\n',maxT8);
assert(minA8>R.alphaPhysicalMin && maxA8<R.alphaPhysicalMax && ...
       maxB8<=B.betaBar+1e-10 && minK8>=B.kappaMin-1e-10 && ...
       maxT8<R.Tmax,'VIB-P8 failed');

Results=struct('VIA',VIA,'R',R,'A',A,'B',B, ...
 'VIB1',struct('maxSO3',maxSO3,'maxAxis',maxAxis,'maxBc',maxBc), ...
 'VIB2',struct('maxComp',maxComp,'maxB',maxB,'maxInertial',maxInertial), ...
 'VIB3',struct('maxTri',maxTri,'maxBeta',maxBeta,'maxAlphaErr',maxAlphaErr, ...
               'maxBetaRatio',maxBetaRatio,'maxAlphaRatio',maxAlphaRatio,'minKappa',minKappa), ...
 'VIB4',struct('minKappa',minK,'maxTdiff',maxTdiff,'maxBoundRatio',maxBoundRatio), ...
 'VIB5',struct('maxBetaDotErr',maxBetaDotErr,'nonInvariantBetaDot',bd0), ...
 'VIB6',struct('alpha',alpha,'beta',beta,'kappa',kappa,'Tphys',Tphys, ...
               'thetaF',thetaF,'thetaR',thetaR,'thetaC',thetaC), ...
 'VIB7',struct('maxWrapC',maxWrapC,'maxStepC',maxStepC,'rawJump',rawJump), ...
 'VIB8',struct('minKappa',minK8,'minAlpha',minA8,'maxAlpha',maxA8, ...
               'maxBeta',maxB8,'maxT',maxT8));
end
