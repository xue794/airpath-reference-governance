function C=cl6_params()
%CL6_PARAMS Complete wind-aware T-2 GTM--GNA 6DOF integration parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% Rev.1 FINAL REGRESSION FREEZE uses the same canonical parameter files
% at every hierarchy level.  No integration-only controller override is
% permitted in this file.

C.revision='Complete Wind-Aware T-2 GTM--GNA 6DOF Closed-Loop Integration -- Final Rev.1';

C.P=t2_gna_params();
C.G=smc_guidance_params();
C.E=wind_estimator_params();
C.R=rg_params();
C.A=viA_params();
C.B=viB_params();
C.F=viC_params();
C.Dp=viD_params();
C.P7=vii_params();

% Full-system simulation tolerances.
C.RelTol=1e-8;
C.AbsTol=1e-10;
C.MaxStep=0.01;
C.sample_dt=0.05;

% Steady windows.
C.tSteadyHold=30.0;
C.tSteadyConstant=20.0;
C.tSteadyVarying=60.0;
C.tSteadyNoisy=60.0;

% Durations.
C.Thold=60.0;
C.Tconstant=40.0;
C.Tvarying=120.0;
C.Tnoisy=120.0;
C.Twrap=60.0;

% Ground-reference geometry.
C.VgHold=C.R.Vb;             % 38.1 m/s
C.VgConstant=C.R.Vb+8.0;     % mean airspeed ~38.1 m/s for w_N=8
C.VgVarying=C.R.Vb+6.0;      % mean airspeed ~38.1 m/s for mean w_N=6
C.wrapSpeed=C.R.Vb;
C.wrapChi0=deg2rad(170);
C.wrapChiRate=0.004;          % rad/s; crosses +pi during 60 s

% Certified simulations initialize exactly on the commanded ground path
% and with the wind estimator initialized to the true wind.  The
% constant_cold case separately tests observer cold start what(0)=0.
C.initialPositionError=zeros(3,1);

% Numerical tolerances used by the automated audit.
C.tol.identity=1e-8;
C.tol.phi=1e-7;
C.tol.qnorm=1e-8;
C.tol.tube=deg2rad(0.70);
C.tol.thetaF=deg2rad(0.30);
C.tol.thetaR=deg2rad(0.40);
end
