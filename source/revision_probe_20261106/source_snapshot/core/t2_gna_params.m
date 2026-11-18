function P = t2_gna_params()
%T2_GNA_PARAMS  T-2 GTM--GNA nominal parameters used in this study.
%
% Rigid-body/geometry values are from Grauer & Morelli (2014), Table 1.
% Aerodynamic coefficients theta(1:45) are the T-2 column of Tables 3--4.
% SI units are used internally. Body axes: x forward, y right, z down.

ft = 0.3048;
ft2 = ft^2;
lbf = 4.4482216152605;
slug = 14.59390293720636;
slugft2 = slug*ft2;

P.revision = 'OpenLoop Rev.1 / WindEstimator Rev.0';
P.g = 9.80665;
P.quat_norm_gain = 0; % pure quaternion kinematics; analytical norm is invariant
P.W_N = 49.6*lbf;
P.m = P.W_N/P.g;
P.Ixx = 1.327*slugft2;
P.Iyy = 4.254*slugft2;
P.Izz = 5.454*slugft2;
P.Ixz = 0.120*slugft2;
% Standard aircraft convention for the product of inertia:
P.I = [P.Ixx 0 -P.Ixz; 0 P.Iyy 0; -P.Ixz 0 P.Izz];
P.S = 5.902*ft2;
P.cbar = 0.915*ft;
P.b = 6.849*ft;
P.xcg_src_ft = 0.25;  % source table entry; not used in this reduced plant

P.theta = [ ...
 +0.019; -0.078; -27.42; +0.293; +3.420; +288.2; -0.040; +1.819; -355.3; -6.563; ...
 -1.003; +0.033; +0.952; -0.009; +0.253; ...
 +0.016; +5.343; +30.78; +0.396; +12.03; +0.506; -36.30; +46.13; ...
 -0.109; -0.366; +0.061; -0.079; +0.021; ...
 +0.182; -1.782; -44.34; -1.785; +374.0; -1748.0; +2.439; +1949.0; -0.038; +0.803; ...
 +0.183; -0.022; -0.405; -0.009; -0.129; +0.184; -0.377];

% Source reference case (Grauer & Morelli Table 5)
P.src.h_m = 1200*ft;
P.src.V_mps = 125*ft;
P.src.alpha_rad = deg2rad(5.0);
P.src.deltae_rad = deg2rad(1.8);
P.src.throttle = 0.12; % reference only; no throttle-to-thrust map is assumed here

% Source modal targets for the GNA case (Table 6)
P.src.mode.spiral = 0.0887;
P.src.mode.phugoid_wn = 0.342;
P.src.mode.phugoid_zeta = 0.0452;
P.src.mode.roll = 5.44;
P.src.mode.dutch_wn = 5.27;
P.src.mode.dutch_zeta = 0.179;
P.src.mode.short_wn = 6.49;
P.src.mode.short_zeta = 0.351;
end
