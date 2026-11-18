function S=simulate_viD_case(mode,T,R,A,F,Dp,Ew)
%SIMULATE_VID_CASE Section V/VI-A + VI-C filter + VI-D ideal rotational plant.

[upMode,distMode]=viD_upstream_mode(mode);
Up=simulate_viA_case(upMode,T,R,A,Ew);

Cmd.t=Up.t;
Cmd.gamma=Up.y(:,2);
Cmd.chi=Up.y(:,3);
Cmd.mu=Up.mu;
Cmd.alpha=Up.alpha;

Cc0=viB_command_attitude(Cmd.gamma(1),Cmd.chi(1),Cmd.mu(1),Cmd.alpha(1));
qr0=viC_dcm_to_quat(Cc0);
omegaR0=zeros(3,1);
Cr0=Cc0;

QR0=rotvec_to_dcm(Dp.initialAxis*Dp.movingTheta0);
C0=Cr0*QR0;
q0=viC_dcm_to_quat(C0);
Arel0=C0.'*Cr0;
omega0=Arel0*omegaR0+Dp.movingEOmega0;

x0=[qr0;omegaR0;q0;omega0];

opts=odeset('RelTol',Dp.RelTol,'AbsTol',Dp.AbsTol,'MaxStep',Dp.MaxStep);
[t,x]=ode45(@(tt,xx)viD_joint_rhs(tt,xx,Up,Cmd,distMode,R,A,F,Dp), ...
            Up.t,x0,opts);

n=numel(t);
thetaF=zeros(n,1); thetaR=zeros(n,1); thetaC=zeros(n,1);
PsiR=zeros(n,1); eR=zeros(n,3); eOmega=zeros(n,3);
tauC=zeros(n,3); dTau=zeros(n,3);
alpha=zeros(n,1); beta=zeros(n,1); kappa=zeros(n,1); Tphys=zeros(n,1);
qNormR=zeros(n,1); qNorm=zeros(n,1);
omegaRDot=zeros(n,3);

for k=1:n
    Cc=viC_interp_command(t(k),Cmd);
    Cr=quat_to_dcm(x(k,1:4).');
    C=quat_to_dcm(x(k,8:11).');
    omegaR=x(k,5:7).';
    omega=x(k,12:14).';

    Qf=Cc.'*Cr;
    eRf=0.5*viC_vee(Qf-Qf.');
    wrDot=-F.kRf*eRf-F.KomegaF*omegaR;
    omegaRDot(k,:)=wrDot.';

    V=Up.y(k,1);
    Bpath=viA_airpath_basis(Up.y(k,2),Up.y(k,3));
    Air=viB_airdata_from_geometry(C,Bpath.t);
    Mom=viD_rate_moment_decomp(V,Air.alpha,R);
    Ctl=viD_virtual_torque(Cr,omegaR,wrDot,C,omega,Mom.D,R,Dp);
    Err=Ctl.Err;

    thetaF(k)=dcm_distance(Cc,Cr);
    thetaR(k)=Err.theta;
    thetaC(k)=dcm_distance(Cc,C);
    PsiR(k)=Err.Psi;
    eR(k,:)=Err.eR.';
    eOmega(k,:)=Err.eOmega.';
    tauC(k,:)=Ctl.tauC.';
    dTau(k,:)=viD_torque_disturbance(t(k),distMode,Dp).';

    alpha(k)=Air.alpha;
    beta(k)=Air.beta;
    kappa(k)=Air.kappa;
    Tphys(k)=Up.Fparallel(k)/Air.kappa;

    qNormR(k)=norm(x(k,1:4));
    qNorm(k)=norm(x(k,8:11));
end

S=struct('t',t,'x',x,'Up',Up,'Cmd',Cmd,'mode',mode, ...
         'thetaF',thetaF,'thetaR',thetaR,'thetaC',thetaC, ...
         'PsiR',PsiR,'eR',eR,'eOmega',eOmega, ...
         'tauC',tauC,'dTau',dTau,'omegaRDot',omegaRDot, ...
         'alpha',alpha,'beta',beta,'kappa',kappa,'Tphys',Tphys, ...
         'qNormR',qNormR,'qNorm',qNorm);
end
