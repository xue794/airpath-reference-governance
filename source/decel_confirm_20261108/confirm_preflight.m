function confirm_preflight(C,O,x0,method,scenario)
calm=O;calm.WindFun=@(t)zeros(6,1);d0_preflight(C,calm,x0,method);
for t=[0 1.5 3 14 20]
 z=O.WindFun(t);a=O.WindFun(t-1e-5);b=O.WindFun(t+1e-5);
 assert(norm((b(1:3)-a(1:3))/2e-5-z(4:6),inf)<1e-8);
 if strcmp(scenario,'C3')&&t>0,assert(norm(z(1:3))>.1);else,assert(all(z==0));end
 Air=d0_air(t,x0(1:13),O.WindFun);
 Cmd=cp6_feedback_command([Air.V;Air.gamma;Air.chi],x0(20:25),C.R,C.A);
 assert(isequal(Air.wI,z(1:3))&&isequal(Air.wDotI,z(4:6)));
 assert(norm(Air.vaI-(Air.vgI-z(1:3)),inf)<1e-10);
 gaps=d0_raw_margins(t,x0,C,O);
 L=rg_N_map(C.R.alphaLo,Air.V,Cmd.X,C.R);U=rg_N_map(C.R.alphaHi,Air.V,Cmd.X,C.R);
 assert(norm(gaps(3:4)-[Cmd.N-L.value;U.value-Cmd.N],inf)<1e-8);
end
end
