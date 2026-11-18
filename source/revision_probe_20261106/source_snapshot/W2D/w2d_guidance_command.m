function G=w2d_guidance_command(t,r,ze,what,mode,C,O)
%W2D_GUIDANCE_COMMAND Frozen Section-IV guidance with only the reference path replaced.
[rd,rdDot,RM]=w2d_reference(t,mode,C,O);
er=r-rd;s=er+C.G.Lambda*ze;sigma=smc_sat(s./C.G.eps);
vc=rdDot-what-C.G.Lambda*er-C.G.Ks*s-C.G.Kc.*sigma;
G.rd=rd;G.rdDot=rdDot;G.er=er;G.s=s;G.sigma=sigma;G.vc=vc;G.ReferenceMeta=RM;
end
