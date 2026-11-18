function Results=run_section_viB_rev0(makeFigures)
%RUN_SECTION_VIB_REV0 Run all VI-B tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_section_viB_rev0();
if makeFigures
    demo_section_viB();
end
end
