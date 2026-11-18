function S = smc_guidance_test_signals(t,mode,E)
%SMC_GUIDANCE_TEST_SIGNALS Exogenous channels for isolated guidance tests.
% Outputs:
%   S.w   true inertial NED wind [m/s]
%   S.wd  true wind derivative [m/s^2]
%   S.eg  governor mismatch v_f-v_c [m/s]
%   S.ev  air-relative velocity realization error v_a-v_f [m/s]
%   S.nw  pseudo-wind measurement error y_w-w [m/s]
%
% The n_w channel is the already-certified estimator interface from WE-P4;
% it is not raw sensor noise.

if nargin < 3
    E = wind_estimator_params();
end

S.eg = zeros(3,1);
S.ev = zeros(3,1);
S.nw = zeros(3,1);

switch lower(mode)
    case 'still'
        S.w = zeros(3,1);
        S.wd = zeros(3,1);

    case 'constant'
        S.w = E.w0;
        S.wd = zeros(3,1);

    case {'varying','mismatch','bound','noisy'}
        [S.w,S.wd] = wind_estimator_test_wind(t,'varying');

        if any(strcmpi(mode,{'mismatch','bound','noisy'}))
            S.eg = [ 0.60*sin(0.17*t);
                    -0.45*cos(0.13*t);
                     0.30*sin(0.11*t+0.20)];
        end

        if any(strcmpi(mode,{'bound','noisy'}))
            S.ev = [0.08*sin(0.31*t);
                    0.06*cos(0.23*t);
                    0.04*sin(0.19*t)];
        end

        if strcmpi(mode,'noisy')
            S.nw = [ 0.28*sin(0.90*t) + 0.08*sin(2.10*t);
                     0.22*cos(0.70*t) - 0.06*sin(1.60*t);
                     0.14*sin(1.10*t+0.40)];
        end

    otherwise
        error('Unknown SMC guidance test mode: %s',mode);
end
end
