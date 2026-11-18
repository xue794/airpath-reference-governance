function Results=test_section_vii_rev0()
%TEST_SECTION_VII_REV0 Section VII tests VII-P1--VII-P8.

VID=test_section_viD_rev0();
R=VID.R;
A=VID.A;
B=VID.B;
F=VID.F;
Dp=VID.Dp;
P7=vii_params();
Ew=wind_estimator_params();

fprintf('\nSection VII Constrained Control Allocation + Actuator Interface Rev.0\n');
fprintf('=====================================================================\n');

%% VII-P1 exact affine moment decomposition and surface-force residual
rng(61,'twister');
maxMoment=0; maxForce=0;
for k=1:500
    V=R.Vmin+(R.Vmax-R.Vmin)*rand;
    alpha=R.alphaPhysicalMin+(R.alphaPhysicalMax-R.alphaPhysicalMin)*rand;
    beta=deg2rad(-0.7+1.4*rand);
    omega=-0.25+0.50*rand(3,1);
    delta=P7.deltaMin+(P7.deltaMax-P7.deltaMin).*rand(3,1);

    va=V*[cos(alpha)*cos(beta);sin(beta);sin(alpha)*cos(beta)];

    [f0,t0]=t2_gna_aero(va,omega,zeros(3,1),eye(3),zeros(3,1),R.h_ref,R.P);
    [f1,t1]=t2_gna_aero(va,omega,delta,eye(3),zeros(3,1),R.h_ref,R.P);

    M=vii_moment_components(V,alpha,beta,R);
    tRec=M.tau0+M.K*omega+M.B*delta;
    fRec=vii_surface_force_residual(V,alpha,delta,R);

    maxMoment=max(maxMoment,norm(t1-tRec));
    maxForce=max(maxForce,norm((f1-f0)-fRec));
end

fprintf('VII-P1 exact T-2 GNA control-affine decomposition:\n');
fprintf('   max aerodynamic moment reconstruction residual = %.3e N m\n',maxMoment);
fprintf('   max surface-force residual reconstruction error = %.3e N\n',maxForce);
assert(maxMoment<1e-11 && maxForce<1e-11,'VII-P1 failed');

%% VII-P2 effectiveness conditioning over certified domain
minSigma=inf; maxCond=0; minAbsDet=inf;
for V=linspace(R.Vmin,R.Vmax,P7.scanNV)
    for alpha=linspace(R.alphaPhysicalMin,R.alphaPhysicalMax,P7.scanNA)
        Be=vii_surface_effectiveness(V,alpha,R);
        sv=svd(Be);
        minSigma=min(minSigma,min(sv));
        maxCond=max(maxCond,max(sv)/min(sv));
        minAbsDet=min(minAbsDet,abs(det(Be)));
    end
end

fprintf('VII-P2 surface-effectiveness conditioning:\n');
fprintf('   grid min sigma_min(B_delta) = %.6f N m/rad\n',minSigma);
fprintf('   grid max condition number   = %.6f\n',maxCond);
fprintf('   grid min |det(B_delta)|     = %.6e\n',minAbsDet);
assert(minSigma>1 && maxCond<10 && minAbsDet>1,'VII-P2 failed');

%% VII-P3 allocator exactness and box-QP optimality
rng(62,'twister');
maxExact=0; maxBoxViolation=0; maxCostGap=0;
for k=1:500
    V=32+12*rand;
    alpha=deg2rad(1+7*rand);
    Be=vii_surface_effectiveness(V,alpha,R);

    dTrue=0.7*P7.deltaMin+1.4*(P7.deltaMax-P7.deltaMin).*rand(3,1);
    dTrue=min(max(dTrue,0.7*P7.deltaMin),0.7*P7.deltaMax);
    req=Be*dTrue;

    O=vii_box_allocator(Be,req,P7.deltaMin,P7.deltaMax,P7.Wtau,P7.allocTol);
    maxExact=max(maxExact,norm(O.residual));
    maxBoxViolation=max(maxBoxViolation,max([P7.deltaMin-O.delta;O.delta-P7.deltaMax;zeros(3,1)]));

    % Saturated stress case: compare against a dense deterministic 31^3 box grid.
    if k<=10
        reqS=(0.8+0.4*rand)*P7.stressTau;
        Os=vii_box_allocator(Be,reqS,P7.deltaMin,P7.deltaMax,P7.Wtau,P7.allocTol);

        grid=linspace(-1,1,31);
        best=inf;
        for ia=1:numel(grid)
            for ie=1:numel(grid)
                for ir=1:numel(grid)
                    dg=[grid(ia)*P7.deltaMax(1);
                        grid(ie)*P7.deltaMax(2);
                        grid(ir)*P7.deltaMax(3)];
                    c=norm(P7.Wtau*(Be*dg-reqS))^2;
                    best=min(best,c);
                end
            end
        end
        maxCostGap=max(maxCostGap,Os.cost-best);
    end
end

fprintf('VII-P3 box-constrained allocator:\n');
fprintf('   max feasible-request residual = %.3e N m\n',maxExact);
fprintf('   max box violation             = %.3e rad\n',maxBoxViolation);
fprintf('   max active-set cost minus dense-grid best = %.3e\n',maxCostGap);
assert(maxExact<1e-10 && maxBoxViolation<1e-12 && maxCostGap<1e-8,'VII-P3 failed');

%% VII-P4 VI-D aggressive virtual-torque allocation, zero servo lag
S4=simulate_viD_case('aggressive',Dp.Taggressive,R,A,F,Dp,Ew);
n=numel(S4.t);
delta4=zeros(n,3); res4=zeros(n,3); active4=zeros(n,3);

for k=1:n
    C=quat_to_dcm(S4.x(k,8:11).');
    omega=S4.x(k,12:14).';
    Pth=viA_airpath_basis(S4.Up.y(k,2),S4.Up.y(k,3));
    Air=viB_airdata_from_geometry(C,Pth.t);

    O=vii_allocate(S4.tauC(k,:).',S4.Up.y(k,1),Air.alpha,Air.beta,omega,R,P7);
    delta4(k,:)=O.deltaCmd.';
    res4(k,:)=O.dAlloc.';
    active4(k,:)=O.active.';
end

maxDef4=max(abs(rad2deg(delta4)),[],1);
maxRes4=max(sqrt(sum(res4.^2,2)));
numActive4=sum(any(active4>0,2));

fprintf('VII-P4 aggressive VI-D torque-request allocation:\n');
fprintf('   max |delta_cmd| [da de dr] = [%.4f %.4f %.4f] deg\n',maxDef4);
fprintf('   max ||d_tau,alloc|| = %.3e N m\n',maxRes4);
fprintf('   active-box samples  = %d / %d\n',numActive4,n);
assert(all(maxDef4<5) && maxRes4<1e-9 && numActive4==0,'VII-P4 failed');

%% VII-P5 standalone actuator position/rate invariance
deltaStep=deg2rad([20;-20;15]);
Tstep=80;
x0=[zeros(3,1);20];

tg=linspace(0,2,2001);
opts=odeset('RelTol',P7.RelTol,'AbsTol',P7.AbsTol,'MaxStep',0.001);
[t5,x5]=ode45(@(tt,xx)vii_step_actuator_rhs(tt,xx,deltaStep,Tstep,P7),tg,x0,opts);

dNum=diff(x5(:,1:3))./diff(t5);
tNum=diff(x5(:,4))./diff(t5);

maxPosViolation=max(max([P7.deltaMin.'-x5(:,1:3);x5(:,1:3)-P7.deltaMax.';zeros(1,3)]));
maxRate=max(max(abs(dNum)));
maxTRate=max(abs(tNum));
Tviol=max([P7.Tmin-min(x5(:,4));max(x5(:,4))-P7.Tmax;0]);

fprintf('VII-P5 servo/propulsor hard constraints:\n');
fprintf('   max surface position violation = %.3e rad\n',maxPosViolation);
fprintf('   max realized surface rate      = %.4f deg/s (limit %.4f)\n', ...
    rad2deg(maxRate),rad2deg(max(P7.deltaRateMax)));
fprintf('   max realized thrust rate       = %.4f N/s (limit %.4f)\n',maxTRate,P7.TRateMax);
fprintf('   thrust-box violation           = %.3e N\n',Tviol);
assert(maxPosViolation<1e-9 && maxRate<=max(P7.deltaRateMax)+1e-5 && ...
       maxTRate<=P7.TRateMax+1e-3 && Tviol<1e-8,'VII-P5 failed');

%% VII-P6 closed rotational loop with physical surface servos
S6=simulate_section_vii_case('aggressive',P7.Taggressive,R,A,F,Dp,P7,Ew);
idx6=S6.t>=P7.tSteady;

maxThetaR6=max(S6.thetaR);
maxTube6=max(S6.thetaF+S6.thetaR);
maxDTau6=max(sqrt(sum(S6.dTau.^2,2)));
lateDTau6=max(sqrt(sum(S6.dTau(idx6,:).^2,2)));
rmsDTau6=sqrt(mean(sum(S6.dTau(idx6,:).^2,2)));
maxDAct6=max(abs(rad2deg(S6.deltaAct)),[],1);
maxDCmd6=max(abs(rad2deg(S6.deltaCmd)),[],1);

fprintf('VII-P6 aggressive closed rotational loop with servos:\n');
fprintf('   max theta_R / total tube = [%.6f %.6f] deg\n', ...
    rad2deg(maxThetaR6),rad2deg(maxTube6));
fprintf('   max / late max / late RMS ||d_tau|| = [%.6f %.6f %.6f] N m\n', ...
    maxDTau6,lateDTau6,rmsDTau6);
fprintf('   max |delta_cmd| [deg] = [%.4f %.4f %.4f]\n',maxDCmd6);
fprintf('   max |delta_act| [deg] = [%.4f %.4f %.4f]\n',maxDAct6);

assert(maxThetaR6<Dp.thetaRMax && maxTube6<B.thetaTube && ...
       lateDTau6<P7.dTauBudget && all(maxDAct6<10),'VII-P6 failed');

%% VII-P7 physical thrust interface + omitted surface-force residual
maxFdelta=max(sqrt(sum(S6.fDelta.^2,2)));
maxFacc=maxFdelta/R.P.m;
maxTerr=max(abs(S6.Tact-S6.Tcmd));
lateTerr=max(abs(S6.Tact(idx6)-S6.Tcmd(idx6)));
Tmin6=min(S6.Tact); Tmax6=max(S6.Tact);

fprintf('VII-P7 physical-force/propulsor interface:\n');
fprintf('   max omitted surface-force residual ||d_F,m,surf|| = %.6f N\n',maxFdelta);
fprintf('   equivalent acceleration magnitude = %.6f m/s^2\n',maxFacc);
fprintf('   max / late max |T_act-T_cmd| = [%.6f %.6f] N\n',maxTerr,lateTerr);
fprintf('   actual thrust range = [%.6f %.6f] N\n',Tmin6,Tmax6);

assert(maxFdelta>0 && Tmin6>=P7.Tmin-1e-8 && Tmax6<=P7.Tmax+1e-8 && ...
       lateTerr<0.25,'VII-P7 failed');

%% VII-P8 saturation stress + topology/wind interface robustness
Be=vii_surface_effectiveness(R.Vb,deg2rad(5),R);
O8=vii_box_allocator(Be,P7.stressTau,P7.deltaMin,P7.deltaMax,P7.Wtau,P7.allocTol);
active8=sum(O8.active);
res8=norm(O8.residual);

Sw=simulate_section_vii_case('wrap',P7.Twrap,R,A,F,Dp,P7,Ew);
Sv=simulate_section_vii_case('varying',P7.Tvarying,R,A,F,Dp,P7,Ew);

rawChi=zeros(numel(Sw.t),1);
for k=1:numel(Sw.t)
    [vc,~]=rg_raw_test_command(Sw.t(k),'wrap',R);
    [cc,~]=rg_raw_command_chart(vc,R);
    rawChi(k)=cc(3);
end
rawJump=max(abs(diff(rawChi)));

maxWrapTube=max(Sw.thetaF+Sw.thetaR);
maxVarTube=max(Sv.thetaF+Sv.thetaR);
lateVarDTau=max(sqrt(sum(Sv.dTau(Sv.t>=P7.tSteady,:).^2,2)));

fprintf('VII-P8 saturation / S1 / wind robustness:\n');
fprintf('   stress allocator active surfaces = %d, residual = %.6f N m\n',active8,res8);
fprintf('   raw course branch jump = %.4f deg\n',rad2deg(rawJump));
fprintf('   wrap max total tube = %.6f deg\n',rad2deg(maxWrapTube));
fprintf('   varying-wind max total tube = %.6f deg\n',rad2deg(maxVarTube));
fprintf('   varying-wind late max ||d_tau|| = %.6f N m\n',lateVarDTau);

assert(active8>=1 && res8>0 && rawJump>deg2rad(300) && ...
       maxWrapTube<B.thetaTube && maxVarTube<B.thetaTube && ...
       lateVarDTau<P7.dTauBudget,'VII-P8 failed');

Results=struct('VID',VID,'R',R,'A',A,'B',B,'F',F,'Dp',Dp,'P7',P7, ...
    'VII1',struct('momentErr',maxMoment,'forceErr',maxForce), ...
    'VII2',struct('minSigma',minSigma,'maxCond',maxCond,'minAbsDet',minAbsDet), ...
    'VII3',struct('maxExact',maxExact,'maxCostGap',maxCostGap), ...
    'VII4',struct('delta',delta4,'residual',res4,'active',active4), ...
    'VII5',struct('t',t5,'x',x5), ...
    'VII6',S6, ...
    'VII7',struct('maxForceResidual',maxFdelta,'maxThrustError',maxTerr), ...
    'VII8',struct('stress',O8,'wrap',Sw,'varying',Sv));
end
