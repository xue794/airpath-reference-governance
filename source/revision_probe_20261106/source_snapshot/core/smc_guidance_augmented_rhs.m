function zDot = smc_guidance_augmented_rhs(t,z,G,E,mode)
%SMC_GUIDANCE_AUGMENTED_RHS Isolated guidance + wind-estimator surrogate.
%
% State:
%   z = [r^I(3); z_e(3); what^I(3)]
%
% Surrogate interfaces:
%   e_g = v_f-v_c
%   e_v = v_a-v_f
%   v_f = v_c+e_g
%   v_a = v_f+e_v
%   rdot = v_a+w
%
% Wind-estimator pseudo measurement:
%   y_w = w+n_w
%
% This deliberately isolates Section IV guidance from Sections V--VII.

r = z(1:3);
ze = z(4:6);
what = z(7:9);

O = smc_guidance_command(t,r,ze,what,G);
S = smc_guidance_test_signals(t,mode,E);

vf = O.vc + S.eg;
va = vf + S.ev;
rDot = va + S.w;

zeDot = O.er - G.Lambda\S.eg;
whatDot = E.Lw*(S.w + S.nw - what);

zDot = [rDot; zeDot; whatDot];
end
