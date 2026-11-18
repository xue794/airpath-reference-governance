function B = w2f_benchmark(R,Bank,O,outdir)
%W2F_BENCHMARK Serial publication timing on the frozen 425-state bank.
% Timings are block averages per governor call. No parallel execution.
N=size(Bank.X,1);
baseOrder=mod((0:N-1)*O.PermutationStride,N)+1;
X=Bank.X.';U0=Bank.U0.';Cbank=Bank.C;uB=Bank.uB;

fprintf('\nWARM-UP -- %d full passes over %d frozen states\n',O.WarmupPasses,N);
for w=1:O.WarmupPasses
    ord=local_shift(baseOrder,(w-1)*O.BlockShift,N);
    local_constraint_batch(X,R,ord,1);
    local_backup_batch(X,R,ord,1);
    local_scalar_selector_batch(U0,Cbank,uB,O,ord,1);
    local_qp_selector_batch(U0,Cbank,uB,O,ord,1);
    local_scalar_whole_batch(X,U0,R,ord,1);
    local_qp_whole_batch(X,U0,R,O,ord,1);
end
pause(O.PauseBeforeTiming_s);

nb=O.NumBlocks;ppb=O.PassesPerBlock;
Block=zeros(nb,10);OrderCell=cell(nb,1);
fprintf('\nMEASURED TIMING -- %d blocks x %d passes/block x %d states\n',nb,ppb,N);
for b=1:nb
    ord=local_shift(baseOrder,(b-1)*O.BlockShift,N);
    tc=local_constraint_batch(X,R,ord,ppb);
    tb=local_backup_batch(X,R,ord,ppb);
    if mod(b,2)==1
        ts=local_scalar_selector_batch(U0,Cbank,uB,O,ord,ppb);
        tq=local_qp_selector_batch(U0,Cbank,uB,O,ord,ppb);
        twq=local_qp_whole_batch(X,U0,R,O,ord,ppb);
        tws=local_scalar_whole_batch(X,U0,R,ord,ppb);
        oc='Selector:S-Q; Whole:Q-S';
    else
        tq=local_qp_selector_batch(U0,Cbank,uB,O,ord,ppb);
        ts=local_scalar_selector_batch(U0,Cbank,uB,O,ord,ppb);
        tws=local_scalar_whole_batch(X,U0,R,ord,ppb);
        twq=local_qp_whole_batch(X,U0,R,O,ord,ppb);
        oc='Selector:Q-S; Whole:S-Q';
    end
    Block(b,:)=[b tc tb ts tq tws twq tq/ts twq/tws (tq-ts)];
    OrderCell{b}=oc;
    fprintf('  block %2d/%2d: selector S=%.6f ms QP=%.6f ms (x%.3f); whole S=%.6f ms QP=%.6f ms (x%.3f)\n', ...
        b,nb,ts,tq,tq/ts,tws,twq,twq/tws);
end

BT=table(Block(:,1),OrderCell,Block(:,2),Block(:,3),Block(:,4),Block(:,5),Block(:,6),Block(:,7),Block(:,8),Block(:,9),Block(:,10), ...
    'VariableNames',{'Block','Order','ConstraintBuild_ms','BackupControl_ms','ScalarSelector_ms','QPSelector_ms', ...
    'ScalarWholeGovernor_ms','QPWholeGovernor_ms','QP_over_ScalarSelector','QP_over_ScalarWhole','QP_minus_ScalarSelector_ms'});
w2f_write_csv(fullfile(outdir,'W2F_timing_blocks.csv'),BT);

names={'ConstraintBuild_ms';'BackupControl_ms';'ScalarSelector_ms';'QPSelector_ms'; ...
       'ScalarWholeGovernor_ms';'QPWholeGovernor_ms';'QP_over_ScalarSelector';'QP_over_ScalarWhole';'QP_minus_ScalarSelector_ms'};
cols={Block(:,2),Block(:,3),Block(:,4),Block(:,5),Block(:,6),Block(:,7),Block(:,8),Block(:,9),Block(:,10)};
M=zeros(numel(names),9);
for k=1:numel(names)
    S=w2f_summary_stats(cols{k});
    M(k,:)=[S.N S.Mean S.Median S.Q05 S.Q25 S.Q75 S.Q95 S.Min S.Max];
end
ST=table(names,M(:,1),M(:,2),M(:,3),M(:,4),M(:,5),M(:,6),M(:,7),M(:,8),M(:,9), ...
    'VariableNames',{'Metric','NBlocks','Mean','Median','Q05','Q25','Q75','Q95','Min','Max'});
w2f_write_csv(fullfile(outdir,'W2F_timing_summary.csv'),ST);

% Publication-facing compact row.
sc=w2f_summary_stats(Block(:,4));qp=w2f_summary_stats(Block(:,5));sw=w2f_summary_stats(Block(:,6));qw=w2f_summary_stats(Block(:,7));
rs=w2f_summary_stats(Block(:,8));rw=w2f_summary_stats(Block(:,9));cc=w2f_summary_stats(Block(:,2));bb=w2f_summary_stats(Block(:,3));
measuredCalls=N*nb*ppb;
Pub=table(N,nb,ppb,measuredCalls,cc.Median,bb.Median,sc.Median,qp.Median,rs.Median,sw.Median,qw.Median,rw.Median, ...
    sc.Q25,sc.Q75,qp.Q25,qp.Q75,sw.Q25,sw.Q75,qw.Q25,qw.Q75, ...
    'VariableNames',{'FrozenStates','Blocks','PassesPerBlock','MeasuredCallsPerTimedItem', ...
    'ConstraintBuildMedian_ms','BackupControlMedian_ms','ScalarSelectorMedian_ms','QPSelectorMedian_ms','QP_over_ScalarSelectorMedian', ...
    'ScalarWholeMedian_ms','QPWholeMedian_ms','QP_over_ScalarWholeMedian', ...
    'ScalarSelectorQ25_ms','ScalarSelectorQ75_ms','QPSelectorQ25_ms','QPSelectorQ75_ms', ...
    'ScalarWholeQ25_ms','ScalarWholeQ75_ms','QPWholeQ25_ms','QPWholeQ75_ms'});
w2f_write_csv(fullfile(outdir,'W2F_publication_timing_table.csv'),Pub);

B=struct('Status','PASS','Blocks',BT,'Summary',ST,'PublicationTable',Pub, ...
    'FrozenStates',N,'NumBlocks',nb,'PassesPerBlock',ppb,'MeasuredCallsPerTimedItem',measuredCalls, ...
    'ProtocolNote','Client-serial, one compute thread, frozen 425-state bank, warmed block averages; no parallel workload.');
end

function ord=local_shift(base,s,N)
s=mod(s,N);ord=circshift(base,[0 s]);
end
function ms=local_constraint_batch(X,R,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);C=rg_constraint_affine(X(:,i),R);sink=sink+C.c(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
function ms=local_backup_batch(X,R,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);u=rg_backup_control(X(:,i),R);sink=sink+u(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
function ms=local_scalar_selector_batch(U0,Cbank,uB,O,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);u=w2f_scalar_selector_core(U0(:,i),Cbank{i},uB(:,i),O.Scale);sink=sink+u(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
function ms=local_qp_selector_batch(U0,Cbank,uB,O,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);u=w2f_qp_selector_core(U0(:,i),Cbank{i},uB(:,i),O);sink=sink+u(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
function ms=local_scalar_whole_batch(X,U0,R,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);G=rg_scalar_governor(X(:,i),U0(:,i),R);sink=sink+G.u(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
function ms=local_qp_whole_batch(X,U0,R,O,ord,passes)
sink=0;z=tic;
for p=1:passes
 for q=1:numel(ord)
  i=ord(q);G=rg_matched_qp(X(:,i),U0(:,i),R,O);sink=sink+G.u(1); %#ok<NASGU>
 end
end
ms=1e3*toc(z)/(passes*numel(ord));
end
