function T=w2c_timing_probe(Cases,O)
%W2C_TIMING_PROBE Development-only, client-serial common-state microbenchmark.
% Not publication timing. Uses Proposed varying/noisy states as a shared bank.
T=struct('Status','SKIPPED','Samples',0,'Repeats',O.TimingRepeats,'ConstraintMedian_ms',NaN, ...
 'ScalarSelectionMedian_ms',NaN,'QPSelectionMedian_ms',NaN,'QP_over_ScalarSelection',NaN, ...
 'ScalarWholeGovernorMedian_ms',NaN,'QPWholeGovernorMedian_ms',NaN,'QP_over_ScalarWhole',NaN, ...
 'Note','DEVELOPMENT_TELEMETRY_NOT_PUBLICATION_BENCHMARK');
if ~O.TimingEnabled,return;end
bankX=[];bankCC=[];modes={'varying','noisy'};
for im=1:numel(modes)
 key=['Proposed_' modes{im}];if ~isfield(Cases,key),continue;end
 S=Cases.(key);n=size(S.Y.xf,1);ns=min(O.TimingSamplesPerStressCase,n);ix=unique(round(linspace(1,n,ns)));
 bankX=[bankX;S.Y.xf(ix,:)];bankCC=[bankCC;S.Y.cc(ix,:)]; %#ok<AGROW>
end
if isempty(bankX),return;end
C0=cl6_params();C0.R.rhoB=O.RhoB;N=size(bankX,1);tc=zeros(N,1);ts=zeros(N,1);tq=zeros(N,1);tws=zeros(N,1);twq=zeros(N,1);
for i=1:N
 x=bankX(i,:).';cc=bankCC(i,:).';u0=rg_nominal_control(x,cc,C0.R);Ccon=rg_constraint_affine(x,C0.R);
 for w=1:O.TimingWarmup,w2c_scalar_from_constraints(x,u0,C0.R,Ccon,O);rg_matched_qp(x,u0,C0.R,O.QP,Ccon);rg_scalar_governor(x,u0,C0.R);rg_matched_qp(x,u0,C0.R,O.QP);end %#ok<NASGU>
 r=O.TimingRepeats;
 z=tic;for q=1:r,rg_constraint_affine(x,C0.R);end;tc(i)=1e3*toc(z)/r;
 z=tic;for q=1:r,w2c_scalar_from_constraints(x,u0,C0.R,Ccon,O);end;ts(i)=1e3*toc(z)/r;
 z=tic;for q=1:r,rg_matched_qp(x,u0,C0.R,O.QP,Ccon);end;tq(i)=1e3*toc(z)/r;
 z=tic;for q=1:r,rg_scalar_governor(x,u0,C0.R);end;tws(i)=1e3*toc(z)/r;
 z=tic;for q=1:r,rg_matched_qp(x,u0,C0.R,O.QP);end;twq(i)=1e3*toc(z)/r;
end
T.Status='COMPLETED';T.Samples=N;T.ConstraintMedian_ms=median(tc);T.ScalarSelectionMedian_ms=median(ts);T.QPSelectionMedian_ms=median(tq);
T.QP_over_ScalarSelection=T.QPSelectionMedian_ms/T.ScalarSelectionMedian_ms;T.ScalarWholeGovernorMedian_ms=median(tws);T.QPWholeGovernorMedian_ms=median(twq);T.QP_over_ScalarWhole=T.QPWholeGovernorMedian_ms/T.ScalarWholeGovernorMedian_ms;
end
