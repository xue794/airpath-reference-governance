function Results = run_section_viA_rev0(makeFigures)
%RUN_SECTION_VIA_REV0 Run all VI-A tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_section_viA_rev0();
if makeFigures
    demo_section_viA('still');
    demo_section_viA('varying');
    demo_section_viA('aggressive');
    demo_section_viA('wrap');
end
end
