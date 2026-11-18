function Results = run_reference_governor_rev0(makeFigures)
%RUN_REFERENCE_GOVERNOR_REV0 Run all Section V tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_reference_governor_rev0();
if makeFigures
    demo_reference_governor('small');
    demo_reference_governor('aggressive');
    demo_reference_governor('wrap');
end
end
