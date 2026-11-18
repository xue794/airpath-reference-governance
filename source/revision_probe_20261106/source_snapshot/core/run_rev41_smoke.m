function Results41Smoke=run_rev41_smoke()
%RUN_REV41_SMOKE Small deterministic preflight before the 500-trial final run.
R=rev41_mc_params();
R.N=20;
R.resume=false;
R.outputDir='rev411_smoke_outputs';
Results41Smoke=run_rev41_monte_carlo(R);
assignin('base','Results41Smoke',Results41Smoke);
end
