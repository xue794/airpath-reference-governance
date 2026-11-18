function E = wind_estimator_params()
%WIND_ESTIMATOR_PARAMS Section IV wind estimator Rev.0 test parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% E.Lw is the FINAL synchronized Rev.1 observer gain used in the paper simulations.

E.revision = 'Section IV Wind Estimator Rev.0';
E.Lw = 2.0*eye(3);            % [1/s], test value only
E.w0 = [8;-4;1];              % [m/s], constant-wind test
E.what0 = [0;0;0];            % [m/s]
E.t_constant = 5.0;           % [s]
E.t_varying = 120.0;          % [s]
E.fit_window = [60 120];      % [s], steady-state sinusoid-fit window

% Deterministic bounded sensor-error test magnitudes.
% n_g bound is on the VECTOR norm. Airdata bounds are componentwise.
E.noise.ng_bar = 0.25;                  % [m/s], ||n_g|| <= ng_bar
E.noise.theta_bar = deg2rad(0.35);      % [rad], attitude principal-angle bound
E.noise.nV_bar = 0.20;                  % [m/s]
E.noise.nalpha_bar = deg2rad(0.20);     % [rad]
E.noise.nbeta_bar  = deg2rad(0.15);     % [rad]
E.noise.nz_bar = norm([E.noise.nV_bar; E.noise.nalpha_bar; E.noise.nbeta_bar]);

% Local airspeed upper bound used only for the theoretical pseudo-wind error bound.
E.noise.Vbar = 45.0;                    % [m/s]
% For T(V,alpha,beta), ||Jacobian||_2 <= max(1,V) on the segment.
E.noise.LT = max(1,E.noise.Vbar + E.noise.nV_bar);
E.noise.nw_bound_generic = E.noise.ng_bar ...
    + 2*E.noise.Vbar*sin(E.noise.theta_bar/2) ...
    + E.noise.LT*E.noise.nz_bar;
% Sharper numerical specialization using componentwise airdata bounds:
% ||Delta T|| <= |nV| + Vbar*(|nalpha|+|nbeta|).
E.noise.nw_bound_componentwise = E.noise.ng_bar ...
    + 2*E.noise.Vbar*sin(E.noise.theta_bar/2) ...
    + E.noise.nV_bar ...
    + E.noise.Vbar*(E.noise.nalpha_bar + E.noise.nbeta_bar);
% Default certificate used by WE-P4.
E.noise.nw_bound = E.noise.nw_bound_componentwise;
end
