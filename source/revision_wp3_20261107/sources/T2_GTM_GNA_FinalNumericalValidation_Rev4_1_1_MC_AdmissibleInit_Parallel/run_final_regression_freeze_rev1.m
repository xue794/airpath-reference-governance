function Results=run_final_regression_freeze_rev1(makeFigures)
%RUN_FINAL_REGRESSION_FREEZE_REV1 Final synchronized acceptance run.
if nargin<1, makeFigures=false; end
Results=test_final_regression_freeze_rev1();
if makeFigures
    demo_complete_6dof_rev0(Results.CL);
end
end
