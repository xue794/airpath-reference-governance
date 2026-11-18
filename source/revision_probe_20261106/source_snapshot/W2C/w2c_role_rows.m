function R=w2c_role_rows(S,O)
%W2C_ROLE_ROWS Per-row limiter/binding/support counts on one trajectory.
Y=S.Y;names=Y.constraintNames;n=size(Y.phi,1);int=sum(Y.intervened);R=cell(numel(names),1);
for j=1:numel(names)
 r=struct('Method',S.Method,'Mode',S.mode,'RowIndex',j,'Row',names{j}, ...
  'Samples',n,'InterventionSamples',int,'AppliedBindingCount',sum(Y.binding(:,j)), ...
  'ScalarLimiterCount',sum(Y.scalarLimiter(:,j)),'QPActiveSetCount',sum(Y.qpActive(:,j)), ...
  'QPSupportCount',sum(Y.qpSupport(:,j)),'MinAppliedPhi',min(Y.phi(:,j)),'MinNominalPhi',min(Y.phi0(:,j)));
 if int>0,r.ScalarLimiterSharePct=100*r.ScalarLimiterCount/int;r.QPSupportSharePct=100*r.QPSupportCount/int;else,r.ScalarLimiterSharePct=0;r.QPSupportSharePct=0;end
 R{j}=r;
end
end
