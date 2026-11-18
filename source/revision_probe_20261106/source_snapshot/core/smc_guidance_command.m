function O = smc_guidance_command(t,r,ze,what,G)
%SMC_GUIDANCE_COMMAND Raw wind-compensated integral-error SMC command.
%
% e_r = r-r_d
% s   = e_r + Lambda z_e
% v_c = rdot_d - what - Lambda e_r - Ks s - Kc sat(s./eps)

[rd,rdDot] = smc_guidance_reference(t);
er = r-rd;
s = er + G.Lambda*ze;
sigma = smc_sat(s./G.eps);

vc = rdDot - what - G.Lambda*er - G.Ks*s - G.Kc.*sigma;

O = struct('rd',rd,'rdDot',rdDot,'er',er,'s',s, ...
           'sigma',sigma,'vc',vc);
end
