function Results=run_complete_6dof_rev1(makeFigures)
%RUN_COMPLETE_6DOF_REV1 Run final synchronized regression and optional figures.
if nargin<1, makeFigures=false; end
Results=test_complete_6dof_rev1();
if makeFigures
    demo_complete_6dof_rev0(Results);
end
end
