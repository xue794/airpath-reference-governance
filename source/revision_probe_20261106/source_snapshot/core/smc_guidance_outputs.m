function Y = smc_guidance_outputs(t,z,G,E,mode)
%SMC_GUIDANCE_OUTPUTS Evaluate all algebraic Section IV guidance signals.
r = z(1:3); ze = z(4:6); what = z(7:9);
O = smc_guidance_command(t,r,ze,what,G);
S = smc_guidance_test_signals(t,mode,E);

vf = O.vc + S.eg;
va = vf + S.ev;
vg = va + S.w;
tildew = S.w-what;

Y = O;
Y.w = S.w;
Y.wd = S.wd;
Y.nw = S.nw;
Y.eg = S.eg;
Y.ev = S.ev;
Y.vf = vf;
Y.va = va;
Y.vg = vg;
Y.tildew = tildew;
end
