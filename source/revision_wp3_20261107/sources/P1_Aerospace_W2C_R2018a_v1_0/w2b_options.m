function O = w2b_options()
%W2B_OPTIONS Frozen W2-B numerical protocol; MATLAB R2018a/base MATLAB.
% These are numerical tolerances, NOT mathematical constraint slack.
O.Version='W2B_R2018a_v1.1';
O.RhoB=0.10;
O.TimeScale_s=1.0;
O.Scale=[2;8*pi/180;12*pi/180]; % m/s^3, rad/s^2, rad/s^2
O.FeasTol=1e-10;             % absolute tolerance after row/variable scaling
O.RelTol=1e-12;
O.DualTol=1e-10;
O.KKTTol=1e-8;
O.RankTol=1e-12;
O.StepTol=1e-12;
O.MaxIterations=100;
O.AdmissionTol=1e-9;
O.ObjectiveTol=1e-8;
O.ReferenceMatchTol=1e-8;
O.Backend='auto';           % auto | enumeration
O.Exhaustive=false;         % reference checks only; not the timed backend
O.FeasibleStart=[];
O.QuadprogChecks=12;        % optional cross-check, not a required toolbox
O.UseParallel=true;         % independent algebra tasks, not timed benchmarks
O.RequestedWorkers=32;      % process workers; user may request 56 explicitly
O.ParallelEnabled=false;    % set only after successful local pool setup
O.ParallelWorkers=0;
O.ParallelReplayChecks=8;   % serial replay of parallel results, same inputs
O.ParallelReplayTol=1e-8;
O.QuadprogAlgorithm='interior-point-convex'; % R2018a-compatible check only
O.QuadprogConstraintTol=1e-10;
O.QuadprogOptimalityTol=1e-10;
O.QuadprogStepTol=1e-12;
O.QuadprogMaxIterations=1000;
O.EnumerationChecks=20;
O.TimingRepeats=2;          % development telemetry, NOT publication timing
end
