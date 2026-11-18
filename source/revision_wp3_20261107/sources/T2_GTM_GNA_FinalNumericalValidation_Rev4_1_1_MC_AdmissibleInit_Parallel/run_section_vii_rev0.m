function Results=run_section_vii_rev0(makeFigures)
%RUN_SECTION_VII_REV0 Run all Section VII tests and optional figures.
if nargin<1, makeFigures=false; end
Results=test_section_vii_rev0();
if makeFigures
    demo_section_vii();
end
end
