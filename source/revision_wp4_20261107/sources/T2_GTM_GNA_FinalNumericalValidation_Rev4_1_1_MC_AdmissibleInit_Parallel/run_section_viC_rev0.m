function Results=run_section_viC_rev0(makeFigures)
%RUN_SECTION_VIC_REV0 Run all Section VI-C tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_section_viC_rev0();
if makeFigures
    demo_section_viC();
end
end
