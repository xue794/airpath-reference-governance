function [x,L]=viii_lift_reduced(t,X,mode,C,P8)
%VIII_LIFT_REDUCED Local lifting X_red -> complete 36-state realization.
%
% This lifting is used only by the Section-VIII numerical differential
% certificate.  It does not alter the controller or simulation model.

if nargin<5 || isempty(P8), P8=viii_params(); end
B=viii_blocks();

z=X(B.idx{1});
s=X(B.idx{2});
wt=X(B.idx{3});
xp=X(B.idx{4});
xG=X(B.idx{5});
xF=X(B.idx{6});
xR=X(B.idx{7});
xA=X(B.idx{8});

K=cl6_case_info(mode,C);
[rd,~]=cl6_reference(t,K.refMode,C);
zw=cl6_wind_profile(t,K.windMode,C);
w=zw(1:3);

% Section-IV variables.
er=s-C.G.Lambda*z;
r=rd+er;
what=w-wt;
Gd=cl6_guidance_command(t,r,z,what,mode,C);
[cc,~]=rg_raw_command_chart(Gd.vc,C.R);

% Governor state from x_G=[eta_c;nu_f].
xf=[cc(1)+xG(1);
    cc(2)+xG(2);
    cc(3)+xG(3);
    xG(4:6)];

% Actual air-path state from xi_p.
V=xf(1)+xp(1);
gamma=xf(2)+xp(2)/P8.Vstar;
chi=xf(3)+xp(3)/P8.Vstar;

if V<=C.R.Vc_min
    error('viii_lift_reduced:V','Lifted airspeed left the local domain.');
end

y=[V;gamma;chi];
Cmd=viA_feedback_command(y,xf,C.R,C.A);
Cc=viB_command_attitude(gamma,chi,Cmd.mu,Cmd.alpha);

% Command-filter attitude from local e_Rf.
Qf=viii_error_vector_to_dcm(xF(1:3));
Cr=Cc*Qf;
omegaR=xF(4:6);

% Actual attitude/rate from local e_R,e_omega.
QR=viii_error_vector_to_dcm(xR(1:3));
CIB=Cr*QR;
Arel=CIB.'*Cr;
omega=Arel*omegaR+xR(4:6);

% Exact ground/air/wind velocity consistency.
Pth=viA_airpath_basis(gamma,chi);
vaI=V*Pth.t;
vgI=vaI+w;
vgB=CIB.'*vgI;

q=viC_dcm_to_quat(CIB);
qr=viC_dcm_to_quat(Cr);

% Exact VI-D/VII commands at the lifted state.
AirGeom=viB_airdata_from_geometry(CIB,Pth.t);
Mom=vii_moment_components(V,AirGeom.alpha,AirGeom.beta,C.R);
eRf=0.5*viC_vee(Qf-Qf.');
omegaRDot=-C.F.kRf*eRf-C.F.KomegaF*omegaR;
Ctl=viD_virtual_torque(Cr,omegaR,omegaRDot,CIB,omega,Mom.D,C.R,C.Dp);
Alloc=cl6_allocate_fast(Ctl.tauC,V,AirGeom.alpha,AirGeom.beta,omega,C);

if AirGeom.kappa<=0
    error('viii_lift_reduced:kappa','Lifted thrust projection is nonpositive.');
end
Tcmd=Cmd.Fparallel/AirGeom.kappa;

delta=Alloc.deltaCmd+xA(1:3);
Tact=Tcmd+xA(4);

if any(delta<C.P7.deltaMin-1e-9) || any(delta>C.P7.deltaMax+1e-9) || ...
        Tact<C.P7.Tmin-1e-9 || Tact>C.P7.Tmax+1e-9
    error('viii_lift_reduced:ActuatorBox', ...
        'Reduced perturbation leaves the physical actuator box.');
end

xpFull=[r;vgB;q;omega];
x=[xpFull;z;what;xf;qr;omegaR;delta;Tact];

L.cc=cc;
L.Cmd=Cmd;
L.Cc=Cc;
L.Cr=Cr;
L.CIB=CIB;
L.deltaCmd=Alloc.deltaCmd;
L.Tcmd=Tcmd;
end
