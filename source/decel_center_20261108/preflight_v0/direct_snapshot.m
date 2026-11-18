function S=direct_snapshot(t,x,mode,C,O)
%W2D_CONTROL_SNAPSHOT Frozen downstream controller with W2-D targeted guidance.
if nargin<6||isempty(O),O=w2d_options();end
xp=x(1:13);r=xp(1:3);ze=x(14:16);what=x(17:19);xf=x(20:25);
qr=x(26:29);omegaR=x(30:32);delta=x(33:35);Tact=x(36);
K=struct('addNoise',false,'mode',mode);
% All W2-D cases use still wind and the same plant model as source hold mode.
Air=d0_air(t,xp,O.WindFun);
Gd=cp6_guidance(t,r,ze,what,C,O);
% The analytic six-component reference was inserted in x upstream of this function.
vfI=rg_airpath_to_velocity(xf(1:3));eg=vfI-Gd.vc;ev=Air.vaI-vfI;
zeDot=Gd.er-C.G.Lambda\eg;
Mwind=wind_estimator_measurement_from_truth(Air.CIB,Air.vgB,Air.wI,t,C.E,K.addNoise);whatDot=C.E.Lw*(Mwind.yw-what);
y=[Air.V;Air.gamma;Air.chi];Cmd=cp6_feedback_command(y,xf,C.R,C.A);Cc=viB_command_attitude(Air.gamma,Air.chi,Cmd.mu,Cmd.alpha);
Cr=quat_to_dcm(qr);Qf=Cc.'*Cr;eRf=0.5*viC_vee(Qf-Qf.');qrDot=quat_kinematics(qr,omegaR,C.F.quatNormGain);omegaRDot=-C.F.kRf*eRf-C.F.KomegaF*omegaR;
Mom=vii_moment_components(Air.V,Air.alpha,Air.beta,C.R);Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot,Air.CIB,Air.omega,Mom.D,C.R,C.Dp);Alloc=cl6_allocate_fast(Ctl.tauC,Air.V,Air.alpha,Air.beta,Air.omega,C);
if Air.kappa<=0,error('W2D:kappa','Body-axis thrust projection lost positivity.');end
Tcmd=Cmd.Fparallel/Air.kappa;Act=vii_actuator_derivative(delta,Tact,Alloc.deltaCmd,Tcmd,C.P7);
S.K=K;S.Air=Air;S.Gd=Gd;S.vfI=vfI;S.eg=eg;S.ev=ev;S.zeDot=zeDot;S.Mwind=Mwind;S.whatDot=whatDot;S.Cmd=Cmd;S.Cc=Cc;S.Cr=Cr;S.Qf=Qf;S.eRf=eRf;S.qrDot=qrDot;S.omegaRDot=omegaRDot;S.Mom=Mom;S.Ctl=Ctl;S.tauC=Ctl.tauC;S.Alloc=Alloc;S.tauReq=Alloc.tauReq;S.deltaCmd=Alloc.deltaCmd;S.Tcmd=Tcmd;S.deltaDot=Act.deltaDot;S.TDot=Act.TDot;S.dAlloc=Alloc.dAlloc;S.dServo=Mom.B*(delta-Alloc.deltaCmd);S.dTauNom=Mom.B*delta-Alloc.tauReq;S.tildew=Air.wI-what;
end
