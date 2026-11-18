function [Rows,Pass]=w2d_witness_reproduction(C,O)
%W2D_WITNESS_REPRODUCTION Reproduce the three predeclared W2-B trim witnesses.
% This is algebraic only; it is not closed-loop evidence.
modes=O.Modes;Rows=cell(numel(modes),1);Pass=true;x=[C.R.Vb;C.R.gammaB;0;0;0;0];
for k=1:numel(modes)
 K=w2d_case_info(modes{k},C,O);cc=K.TargetAirPath;u0=rg_nominal_control(x,cc,C.R);Gp=w2d_governor(x,u0,C.R,'Proposed',O,[]);Gq=w2d_governor(x,u0,C.R,'MatchedQP',O,[]);
 lm=min(Gp.lambdaJ);lt=O.LimiterTol*(1+abs(lm));lim=find(abs(Gp.lambdaJ-lm)<=lt);sup=Gq.Solver.SupportRows(:);idx=K.ExpectedRowIndex;
 hitP=any(lim==idx);hitQ=any(sup==idx);row=struct('Mode',modes{k},'ExpectedRowIndex',idx,'ExpectedRow',K.ExpectedRow,'RawV',cc(1),'RawGamma_deg',rad2deg(cc(2)),'RawChi_deg',rad2deg(cc(3)),'VB_over_rho',0,'ProposedLambda',Gp.lambda,'ProposedObjective',Gp.Objective,'QPObjective',Gq.Objective,'ProposedExpectedLimiter',hitP,'QPExpectedSupport',hitQ,'QPPath',Gq.Solver.Path,'MinAppliedPhiQP',min(Gq.phi));Rows{k}=row;Pass=Pass&&hitP&&hitQ&&Gq.Success;
end
end
