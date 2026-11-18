function Results=run_section_viD_rev0(makeFigures)
%RUN_SECTION_VID_REV0 Run all VI-D tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_section_viD_rev0();
if makeFigures
    demo_section_viD();
end
end
