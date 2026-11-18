function confirm_preflight(C,O,x0,method,scenario)
calm=O;calm.WindFun=@(t)zeros(6,1);d0_preflight(C,calm,x0,method);
for t=[0 1.5 3 14 20]
 z=O.WindFun(t);a=O.WindFun(t-1e-5);b=O.WindFun(t+1e-5);
 assert(norm((b(1:3)-a(1:3))/2e-5-z(4:6),inf)<1e-8);
 if strcmp(scenario,'C3')&&t>0,assert(norm(z(1:3))>.1);else,assert(all(z==0));end
 S=d0_snapshot(t,x0,'turn',C,method,O);
 assert(isequal(S.Air.wI,z(1:3))&&isequal(S.Air.wDotI,z(4:6)));
 assert(norm(S.Air.vaI-(S.Air.vgI-z(1:3)),inf)<1e-10);
 gaps=d0_raw_margins(t,x0,C,O);
 L=rg_N_map(C.R.alphaLo,S.Air.V,S.Cmd.X,C.R);U=rg_N_map(C.R.alphaHi,S.Air.V,S.Cmd.X,C.R);
 assert(norm(gaps(3:4)-[S.Cmd.N-L.value;U.value-S.Cmd.N],inf)<1e-8);
end
end
