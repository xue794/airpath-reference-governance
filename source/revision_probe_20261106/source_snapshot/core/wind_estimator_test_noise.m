function N = wind_estimator_test_noise(t,E)
%WIND_ESTIMATOR_TEST_NOISE Smooth deterministic bounded measurement errors.
% The construction is deterministic so adaptive ODE solvers remain repeatable.
%
% Rev.4.1: when E.MC is present, the waveform parameters are randomized once
% per Monte-Carlo run and then remain deterministic during that ODE solve.
% The randomized waveform never exceeds the original componentwise bounds.

if isfield(E,'MC')
    M=E.MC;
    sc=M.noiseScale;

    % Ground-velocity error, vector norm <= sc*ng_bar <= ng_bar.
    d=[sin(0.71*M.freqScale(1)*t+M.phase(1)); ...
       cos(0.43*M.freqScale(2)*t+M.phase(2)); ...
       0.6*sin(0.31*M.freqScale(3)*t+M.phase(3))];
    if norm(d)>0, d=d/norm(d); end
    ampg=E.noise.ng_bar*sc*(0.65+0.35*sin(0.17*M.freqScale(4)*t+M.phase(4))^2);
    N.ng=ampg*d;

    % Attitude measurement error, bounded rotation-vector magnitude.
    a=[cos(0.29*M.freqScale(5)*t+M.phase(5)); ...
       sin(0.37*M.freqScale(6)*t+M.phase(6)); ...
       0.5*cos(0.19*M.freqScale(7)*t+M.phase(7))];
    if norm(a)>0, a=a/norm(a); end
    th=E.noise.theta_bar*sc*(0.55+0.45*sin(0.23*M.freqScale(8)*t+M.phase(8))^2);
    N.dtheta=th*a;

    % Airdata errors stay within the original componentwise magnitudes.
    N.nz=[sc*E.noise.nV_bar*sin(0.83*M.freqScale(9)*t+M.phase(9)); ...
          sc*E.noise.nalpha_bar*cos(0.59*M.freqScale(10)*t+M.phase(10)); ...
          sc*E.noise.nbeta_bar*sin(0.47*M.freqScale(11)*t+M.phase(11))];
    return;
end

% Original deterministic Rev.4 waveform.
d=[sin(0.71*t+0.2); cos(0.43*t-0.1); 0.6*sin(0.31*t+0.7)];
if norm(d)>0, d=d/norm(d); end
ampg=E.noise.ng_bar*(0.65+0.35*sin(0.17*t)^2);
N.ng=ampg*d;

a=[cos(0.29*t); sin(0.37*t+0.4); 0.5*cos(0.19*t-0.2)];
if norm(a)>0, a=a/norm(a); end
th=E.noise.theta_bar*(0.55+0.45*sin(0.23*t+0.1)^2);
N.dtheta=th*a;

N.nz=[E.noise.nV_bar*sin(0.83*t); ...
      E.noise.nalpha_bar*cos(0.59*t+0.3); ...
      E.noise.nbeta_bar*sin(0.47*t-0.5)];
end
