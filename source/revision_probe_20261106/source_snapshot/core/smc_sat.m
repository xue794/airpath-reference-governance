function y = smc_sat(x)
%SMC_SAT Componentwise unit saturation.
y = min(max(x,-1),1);
end
