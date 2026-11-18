function [x0,I0]=cl6_initial_state(mode,C)
%CL6_INITIAL_STATE Consistent local initial state for the full 6DOF loop.
%
% Certified cases begin on the desired ground path with actual air-path
% state equal to the governed reference.  The constant_cold case differs
% only by setting what(0)=0, intentionally exercising observer startup.

K=cl6_case_info(mode,C);
[rd,rdDot]=cl6_reference(0,K.refMode,C);
zW=cl6_wind_profile(0,K.windMode,C);
w0=zW(1:3);

% Raw air-relative command implied by the desired ground velocity and wind.
vaIraw0=rdDot-w0;
Vraw0=norm(vaIraw0);
Vhraw0=hypot(vaIraw0(1),vaIraw0(2));
gammaRaw0=atan2(-vaIraw0(3),Vhraw0);
chiRaw0=atan2(vaIraw0(2),vaIraw0(1));

% Rev.4.1.1 Monte-Carlo admissible initialization.
%
% A constrained reference governor requires its internal state to start in
% the declared backup-feasible set.  The Rev.4.1 pilot instead initialized
% x_f at the random raw air-relative command; 287/500 samples were therefore
% outside the backup set at t=0.  Rev.4.1.1 keeps the SAME uncertainty draws
% and controller gains but starts x_f at the strict backup trim.  The raw
% command remains unchanged and is handled by the governor after t=0.
if isfield(C,'MC') && isfield(C.MC,'initStrategy') && ...
        strcmp(C.MC.initStrategy,'backup_feasible')
    V0=C.R.Vb;
    gamma0=C.R.gammaB;
    chi0=chiRaw0;

    vaI0=V0*[cos(gamma0)*cos(chi0); ...
             cos(gamma0)*sin(chi0); ...
            -sin(gamma0)];
else
    vaI0=vaIraw0;
    V0=Vraw0;
    gamma0=gammaRaw0;
    chi0=chiRaw0;
end

xf0=[V0;gamma0;chi0;0;0;0];
y0=xf0(1:3);

Cmd0=viA_feedback_command(y0,xf0,C.R,C.A);
Cc0=viB_command_attitude(gamma0,chi0,Cmd0.mu,Cmd0.alpha);

% Start both actual and filtered attitudes exactly at Cc0.
CIB0=Cc0;
Cr0=Cc0;
q0=viC_dcm_to_quat(CIB0);
qr0=viC_dcm_to_quat(Cr0);

omega0=zeros(3,1);
omegaR0=zeros(3,1);
omegaRDot0=zeros(3,1);

% Exact air/ground velocity consistency.
vaB0=CIB0.'*vaI0;
vgB0=vaB0+CIB0.'*w0;

% Initial surfaces exactly realize the zero-error rotational virtual port.
Mom0=vii_moment_components(V0,Cmd0.alpha,0,C.R);
Ctl0=viD_virtual_torque(Cr0,omegaR0,omegaRDot0, ...
    CIB0,omega0,Mom0.D,C.R,C.Dp);
Alloc0=cl6_allocate_fast(Ctl0.tauC,V0,Cmd0.alpha,0,omega0,C);

T0=Cmd0.Fparallel/cos(Cmd0.alpha);
T0=min(max(T0,C.P7.Tmin),C.P7.Tmax);

r0=rd+C.initialPositionError;
ze0=zeros(3,1);

if K.coldWindStart
    what0=zeros(3,1);
else
    what0=w0;
end

xp0=[r0;vgB0;q0;omega0];
x0=[xp0;ze0;what0;xf0;qr0;omegaR0;Alloc0.deltaCmd;T0];

% Rev.4.1 Monte-Carlo initial-condition perturbations.
% Existing deterministic modes do not enter this branch.
if isfield(C,'MC') && isfield(C.MC,'init')
    J=C.MC.init;

    x0(1:3)=x0(1:3)+J.position(:);
    x0(4:6)=x0(4:6)+J.vgB(:);

    CIBmc=CIB0*rotvec_to_dcm(J.attitudeRotvec(:));
    x0(7:10)=viC_dcm_to_quat(CIBmc);

    x0(11:13)=x0(11:13)+J.omega(:);
    x0(17:19)=w0+J.windEstimateError(:);
end

I0.rd=rd;
I0.rdDot=rdDot;
I0.w0=w0;
I0.vaI0=vaI0;
I0.vaIraw0=vaIraw0;
I0.Vraw0=Vraw0;
I0.gammaRaw0=gammaRaw0;
I0.chiRaw0=chiRaw0;
I0.V0=V0;
I0.gamma0=gamma0;
I0.chi0=chi0;
I0.Cmd0=Cmd0;
I0.Cc0=Cc0;
I0.delta0=Alloc0.deltaCmd;
I0.T0=T0;
I0.what0=what0;
end
