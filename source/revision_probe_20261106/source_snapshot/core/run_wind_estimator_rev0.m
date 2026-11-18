function Results = run_wind_estimator_rev0(makePlots)
%RUN_WIND_ESTIMATOR_REV0 Run WE-P1--P5 and optional estimator demo.
if nargin<1, makePlots=false; end
Results=test_wind_estimator_rev0();
if makePlots
    demo_wind_estimator(false);
    demo_wind_estimator(true);
end
end
