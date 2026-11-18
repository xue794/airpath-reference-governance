function G=cp6_guidance(t,r,ze,what,C,O)
[rd,rdDot]=cp6_reference(t,C,O);
er=r-rd;s=er+C.G.Lambda*ze;sigma=smc_sat(s./C.G.eps);
vc=rdDot-what-C.G.Lambda*er-C.G.Ks*s-C.G.Kc.*sigma;
G=struct('rd',rd,'rdDot',rdDot,'er',er,'s',s,'sigma',sigma,'vc',vc,'ReferenceMeta',[]);
end
