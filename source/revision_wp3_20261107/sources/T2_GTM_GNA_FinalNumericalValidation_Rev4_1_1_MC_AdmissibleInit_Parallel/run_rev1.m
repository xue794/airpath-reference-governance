function Results = run_rev1(runDemos)
%RUN_REV1 Run Rev.1 automated tests and optionally the demonstration plots.
if nargin<1, runDemos=true; end
clc;
Results=test_open_loop_rev1();
if runDemos
    demo_open_loop();
    demo_quat_equivalence();
    demo_wind_kinematics();
    demo_time_varying_wind();
end
end
