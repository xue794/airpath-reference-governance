function R = rg_params()
%RG_PARAMS Section V Certified Constrained Air-Path Reference Governor Rev.0.
% MATLAB R2018a compatible; base MATLAB only.
%
% FINAL Rev.1 synchronized local-certificate / simulation values.

P = t2_gna_params();

R.revision = 'Section V Certified Constrained Air-Path Reference Governor Rev.0';
R.P = P;
R.g = P.g;
R.h_ref = P.src.h_m;
R.rho = t2_isa_density(R.h_ref);

% Raw-command chart regularity.
R.Vc_min = 5.0;             % m/s
R.Vhc_min = 2.0;            % m/s
R.course_chart_margin = deg2rad(20);

% Second-order reference generator.
R.Zeta = diag([0.90 0.90 0.90]);
R.Omega = diag([1.00 1.20 1.00]);

% Reference-state constraints.
R.Vmin = 30.0;              % m/s
R.Vmax = 46.0;              % m/s
R.gammaMax = deg2rad(10);   % rad
R.nuMax = [2.0; deg2rad(8); deg2rad(12)];

% Normal-force / bank / AoA / thrust feasibility.
R.aGammaMin = 4.0;          % m/s^2
R.muMax = deg2rad(35);
R.alphaPhysicalMin = deg2rad(-2);
R.alphaPhysicalMax = deg2rad(10);
R.alphaLo = deg2rad(0);     % certified inversion interval
R.alphaHi = deg2rad(9);
R.DeltaN = 10.0;            % N
R.Tmin = 0.0;               % N
R.Tmax = 60.0;              % N
R.DeltaT = 3.0;             % N

% CBF/HOCBF gains.
R.k1 = 1.20;
R.k2 = 1.20;
R.kRate = 1.50;
R.kAGamma = 1.00;
R.kMu = 1.00;
R.kN = 1.00;
R.kT = 1.00;
R.kB = 0.40;

% Strict backup trim.
R.Vb = P.src.V_mps;         % 38.1 m/s
R.gammaB = 0.0;
R.kVp = 0.60;
R.kVd = 1.40;
R.kGammaP = 1.20;
R.kGammaD = 1.80;
R.kChiD = 1.50;

R.AB = [0 0 1 0 0;
        0 0 0 1 0;
       -R.kVp 0 -R.kVd 0 0;
        0 -R.kGammaP 0 -R.kGammaD 0;
        0 0 0 0 -R.kChiD];

% Q_B is deliberately dimensionally weighted.  The theory permits any Q_B>0.
sc = [5.0; deg2rad(5); 1.0; deg2rad(5); deg2rad(8)];
R.QB = diag(1./(sc.^2));
R.PB = rg_solve_lyap(R.AB,R.QB);
R.rhoB = 0.20;

% Module-test simulation settings.
R.RelTol = 1e-9;
R.AbsTol = 1e-11;
R.Tsmall = 60;
R.Taggressive = 120;
R.Twrap = 50;
R.sample_dt = 0.05;
end
