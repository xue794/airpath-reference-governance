function d0_preflight(C,O,x0,method)
% One fixed nonmission wind fixture exercises the plant and observer/event wiring.
t=2;xf=x0(20:25);base=d0_snapshot(t,x0,'turn',C,method,O);
assert(isequal(d0_raw_margins(t,x0,C,O),cp6_raw_margins(t,x0,C)));
O.WindFun=@(s)[.2*s;-.1*s;.03*s;.2;-.1;.03];
test=d0_snapshot(t,x0,'turn',C,method,O);
U=struct('delta',x0(33:35),'T',x0(36),'wind_fun',O.WindFun);
[w,wd]=t2_wind_eval(t,U); assert(isequal(test.Air.wI,w)&&isequal(test.Air.wDotI,wd));
dp=t2_gna_6dof_quat(t,x0(1:13),U,C.P);
U.wind_fun=@(s)zeros(6,1);dp0=t2_gna_6dof_quat(t,x0(1:13),U,C.P);
assert(norm(dp-dp0)>1e-5 && norm(test.whatDot-base.whatDot)>1e-3);
gaps=d0_raw_margins(t,x0,C,O);
assert(norm(gaps(3:4)-[test.Cmd.N-rg_force(C.R.alphaLo,test,C);rg_force(C.R.alphaHi,test,C)-test.Cmd.N],inf)<1e-8);
assert(norm(test.Air.vaI-(test.Air.vgI-w),inf)<1e-10);
end
function v=rg_force(alpha,S,C)
M=rg_N_map(alpha,S.Air.V,S.Cmd.X,C.R);v=M.value;
end
