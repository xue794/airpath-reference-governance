function Results = run_smc_guidance_rev0(makeFigures)
%RUN_SMC_GUIDANCE_REV0 Run all Section IV guidance tests and optional demos.
if nargin<1, makeFigures=false; end
Results=test_smc_guidance_rev0();
if makeFigures
    demo_smc_guidance('varying');
    demo_smc_guidance('bound');
    demo_smc_guidance('noisy');
end
end
