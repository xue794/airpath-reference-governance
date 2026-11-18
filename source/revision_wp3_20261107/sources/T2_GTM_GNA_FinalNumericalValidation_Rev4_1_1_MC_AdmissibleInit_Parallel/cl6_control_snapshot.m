function S=cl6_control_snapshot(t,x,mode,C)
%CL6_CONTROL_SNAPSHOT All algebraic controller/interface quantities.
%
% State ordering:
%  1:13  T-2 quaternion 6DOF plant
% 14:16  z_e
% 17:19  wind estimate what^I
% 20:25  second-order governed air-path reference x_f
% 26:29  command-filter quaternion q_r
% 30:32  command-filter angular rate omega_r
% 33:35  actual surface deflections
% 36     actual thrust

xp=x(1:13);
r=xp(1:3);
ze=x(14:16);
what=x(17:19);
xf=x(20:25);
qr=x(26:29);
omegaR=x(30:32);
delta=x(33:35);
Tact=x(36);

K=cl6_case_info(mode,C);
Air=cl6_airpath_state(t,xp,mode,C);

% Section IV SMC guidance.
Gd=cl6_guidance_command(t,r,ze,what,mode,C);

% Section V certified second-order reference governor.
[cc,chart]=rg_raw_command_chart(Gd.vc,C.R);
u0=rg_nominal_control(xf,cc,C.R);
Gov=rg_scalar_governor(xf,u0,C.R);
xfDot=[xf(4:6);Gov.u];

vfI=rg_airpath_to_velocity(xf(1:3));
eg=vfI-Gd.vc;
ev=Air.vaI-vfI;

zeDot=Gd.er-C.G.Lambda\eg;

% Section IV wind estimator, now driven by actual 6DOF truth measurements.
Mwind=wind_estimator_measurement_from_truth( ...
    Air.CIB,Air.vgB,Air.wI,t,C.E,K.addNoise);
whatDot=C.E.Lw*(Mwind.yw-what);

% Section VI-A feedback force/path command from the ACTUAL air-path state.
y=[Air.V;Air.gamma;Air.chi];
Cmd=viA_feedback_command(y,xf,C.R,C.A);

% Section VI-B/VI-C command attitude and dynamic command filter.
Cc=viB_command_attitude(Air.gamma,Air.chi,Cmd.mu,Cmd.alpha);
Cr=quat_to_dcm(qr);

Qf=Cc.'*Cr;
eRf=0.5*viC_vee(Qf-Qf.');

qrDot=quat_kinematics(qr,omegaR,C.F.quatNormGain);
omegaRDot=-C.F.kRf*eRf-C.F.KomegaF*omegaR;

% Section VI-D rotational virtual torque, using actual C,omega.
Mom=vii_moment_components(Air.V,Air.alpha,Air.beta,C.R);
Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot, ...
    Air.CIB,Air.omega,Mom.D,C.R,C.Dp);

% Section VII allocation and actuator interface.
Alloc=cl6_allocate_fast(Ctl.tauC,Air.V,Air.alpha,Air.beta,Air.omega,C);

if Air.kappa<=0
    error('cl6_control_snapshot:kappa','Body-axis thrust projection lost positivity.');
end
Tcmd=Cmd.Fparallel/Air.kappa;

Act=vii_actuator_derivative(delta,Tact,Alloc.deltaCmd,Tcmd,C.P7);

S.K=K;
S.Air=Air;
S.Gd=Gd;
S.cc=cc;
S.chart=chart;
S.Gov=Gov;
S.xfDot=xfDot;
S.vfI=vfI;
S.eg=eg;
S.ev=ev;
S.zeDot=zeDot;
S.Mwind=Mwind;
S.whatDot=whatDot;
S.Cmd=Cmd;
S.Cc=Cc;
S.Cr=Cr;
S.Qf=Qf;
S.eRf=eRf;
S.qrDot=qrDot;
S.omegaRDot=omegaRDot;
S.Mom=Mom;
S.Ctl=Ctl;
S.tauC=Ctl.tauC;
S.Alloc=Alloc;
S.tauReq=Alloc.tauReq;
S.deltaCmd=Alloc.deltaCmd;
S.Tcmd=Tcmd;
S.deltaDot=Act.deltaDot;
S.TDot=Act.TDot;
S.dAlloc=Alloc.dAlloc;
S.dServo=Mom.B*(delta-Alloc.deltaCmd);
S.dTauNom=Mom.B*delta-Alloc.tauReq;
S.tildew=Air.wI-what;
end
