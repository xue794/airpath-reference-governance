function S=simulate_section_vii_case(mode,T,R,A,F,Dp,P7,Ew)
%SIMULATE_SECTION_VII_CASE Dynamic physical allocation + actuator interface.

if strcmpi(mode,'aggressive')
    upMode='aggressive';
elseif strcmpi(mode,'varying')
    upMode='varying';
elseif strcmpi(mode,'wrap')
    upMode='wrap';
else
    error('Unknown Section VII mode: %s',mode);
end

Up=simulate_viA_case(upMode,T,R,A,Ew);
Cmd.t=Up.t;
Cmd.gamma=Up.y(:,2);
Cmd.chi=Up.y(:,3);
Cmd.mu=Up.mu;
Cmd.alpha=Up.alpha;

% Initial VI-C / VI-D states.
Cc0=viB_command_attitude(Cmd.gamma(1),Cmd.chi(1),Cmd.mu(1),Cmd.alpha(1));
qr0=viC_dcm_to_quat(Cc0);
omegaR0=zeros(3,1);
Cr0=Cc0;

QR0=rotvec_to_dcm(Dp.initialAxis*Dp.movingTheta0);
C0=Cr0*QR0;
q0=viC_dcm_to_quat(C0);
omega0=C0.'*Cr0*omegaR0+Dp.movingEOmega0;

% Initialize surfaces/thrust at the exact current allocation to avoid an
% artificial "servo starts at zero" trim impulse.
Pth0=viA_airpath_basis(Up.y(1,2),Up.y(1,3));
Air0=viB_airdata_from_geometry(C0,Pth0.t);
M0=vii_moment_components(Up.y(1,1),Air0.alpha,Air0.beta,R);
eRf0=0.5*viC_vee(Cc0.'*Cr0-(Cc0.'*Cr0).');
omegaRDot0=-F.kRf*eRf0-F.KomegaF*omegaR0;
Ctl0=viD_virtual_torque(Cr0,omegaR0,omegaRDot0,C0,omega0,M0.D,R,Dp);
Alloc0=vii_allocate(Ctl0.tauC,Up.y(1,1),Air0.alpha,Air0.beta,omega0,R,P7);

Tcmd0=Up.Fparallel(1)/Air0.kappa;
T0=min(max(Tcmd0,P7.Tmin),P7.Tmax);

x0=[qr0;omegaR0;q0;omega0;Alloc0.deltaCmd;T0];

opts=odeset('RelTol',P7.RelTol,'AbsTol',P7.AbsTol,'MaxStep',P7.MaxStep);
[t,x]=ode45(@(tt,xx)vii_joint_rhs(tt,xx,Up,Cmd,R,A,F,Dp,P7), ...
            Up.t,x0,opts);

n=numel(t);
thetaF=zeros(n,1); thetaR=zeros(n,1);
alpha=zeros(n,1); beta=zeros(n,1); kappa=zeros(n,1);
deltaCmd=zeros(n,3); deltaAct=x(:,15:17);
Tcmd=zeros(n,1); Tact=x(:,18);
dAlloc=zeros(n,3); dServo=zeros(n,3); dTau=zeros(n,3);
tauReq=zeros(n,3); tauC=zeros(n,3); tauPhys=zeros(n,3);
fDelta=zeros(n,3); active=zeros(n,3);

for k=1:n
    Cc=viC_interp_command(t(k),Cmd);
    Cr=quat_to_dcm(x(k,1:4).');
    C=quat_to_dcm(x(k,8:11).');
    omegaR=x(k,5:7).';
    omega=x(k,12:14).';

    Qf=Cc.'*Cr;
    eRf=0.5*viC_vee(Qf-Qf.');
    wrDot=-F.kRf*eRf-F.KomegaF*omegaR;

    Pth=viA_airpath_basis(Up.y(k,2),Up.y(k,3));
    Air=viB_airdata_from_geometry(C,Pth.t);
    Mom=vii_moment_components(Up.y(k,1),Air.alpha,Air.beta,R);
    Ctl=viD_virtual_torque(Cr,omegaR,wrDot,C,omega,Mom.D,R,Dp);
    Alloc=vii_allocate(Ctl.tauC,Up.y(k,1),Air.alpha,Air.beta,omega,R,P7);

    tc=Up.Fparallel(k)/Air.kappa;

    thetaF(k)=dcm_distance(Cc,Cr);
    thetaR(k)=dcm_distance(Cr,C);
    alpha(k)=Air.alpha;
    beta(k)=Air.beta;
    kappa(k)=Air.kappa;

    deltaCmd(k,:)=Alloc.deltaCmd.';
    dAlloc(k,:)=Alloc.dAlloc.';
    dServo(k,:)=(Mom.B*(deltaAct(k,:).'-Alloc.deltaCmd)).';
    dTau(k,:)=(Mom.B*deltaAct(k,:).'-Alloc.tauReq).';

    tauReq(k,:)=Alloc.tauReq.';
    tauC(k,:)=Ctl.tauC.';
    tauPhys(k,:)=(Mom.tau0+Mom.K*omega+Mom.B*deltaAct(k,:).').';

    fDelta(k,:)=vii_surface_force_residual(Up.y(k,1),Air.alpha,deltaAct(k,:).',R).';
    active(k,:)=Alloc.active.';

    Tcmd(k)=tc;
end

S=struct('t',t,'x',x,'Up',Up,'Cmd',Cmd, ...
    'thetaF',thetaF,'thetaR',thetaR,'alpha',alpha,'beta',beta,'kappa',kappa, ...
    'deltaCmd',deltaCmd,'deltaAct',deltaAct,'Tcmd',Tcmd,'Tact',Tact, ...
    'dAlloc',dAlloc,'dServo',dServo,'dTau',dTau, ...
    'tauReq',tauReq,'tauC',tauC,'tauPhys',tauPhys, ...
    'fDelta',fDelta,'active',active,'mode',mode);
end
