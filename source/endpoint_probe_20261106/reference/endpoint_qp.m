function [u,D,C]=endpoint_qp(t,x,S,R,O,method,allowExtension)
if nargin<7,allowExtension=false;end
RR=R;if strcmp(method,'full17_N0'),RR.DeltaN=0;end
if allowExtension,C=endpoint_affine_extended(x,RR);else,C=rg_constraint_affine(x,RR);end
if strcmp(method,'thrust15'),ix=[1:13 16:17];else,ix=1:17;end
u0=rg_nominal_control(x,endpoint_profile(t,S,RR),RR);Q=O;Q.FeasibleStart=rg_backup_control(x,RR);
[u,D]=w2b_qp3(u0,-C.G(ix,:),C.c(ix),Q.Scale,Q);
if ~D.Success,error('EndpointProbe:QPNoVerifiedSolution','%s',D.Status);end
D.rowIndices=ix;D.supportOriginal=ix(D.SupportRows);
end
