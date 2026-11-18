function O=cl6_guidance_command(t,r,ze,what,mode,C)
%CL6_GUIDANCE_COMMAND Full-system Section-IV SMC guidance command.

K=cl6_case_info(mode,C);
[rd,rdDot]=cl6_reference(t,K.refMode,C);

er=r-rd;
s=er+C.G.Lambda*ze;
sigma=smc_sat(s./C.G.eps);

vc=rdDot-what-C.G.Lambda*er-C.G.Ks*s-C.G.Kc.*sigma;

O.rd=rd;
O.rdDot=rdDot;
O.er=er;
O.s=s;
O.sigma=sigma;
O.vc=vc;
end
