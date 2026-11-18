function audit_interval_scope()
% Counterfactual root existence under code's wider alpha diagnostic interval.
% This does not change any running controller or certify the physical model.
here=fileparts(mfilename('fullpath'));
addpath('D:/Research/paper1_Submission/revision_probe_20261106/source_snapshot/core');
R=rg_params();R.rhoB=.10;
X=[31 0 0 0 0 8;31 0 0 0 4 0;31 0 0 0 0 0];
X(:,[2 3 5 6])=X(:,[2 3 5 6])*pi/180;
labels={'turn8degps','pullup4degps','level31'}; A=struct([]);
for k=1:size(X,1)
 x=X(k,:)';F=rg_feedforward(x,R);
 E=R;E.alphaLo=R.alphaPhysicalMin;E.alphaHi=R.alphaPhysicalMax;
 FE=rg_feedforward(x,E);
 q=struct('name',labels{k},'state',x','original_interval_deg',[0 9], ...
  'original_bracket',F.rootInfo.bracket,'original_alpha_deg',F.alpha*180/pi, ...
  'original_raw_upper_N',F.hNplus+R.DeltaN, ...
  'extended_interval_deg',[E.alphaLo E.alphaHi]*180/pi, ...
  'extended_bracket',FE.rootInfo.bracket,'extended_alpha_deg',FE.alpha*180/pi, ...
  'extended_Tff_N',FE.T,'extended_raw_upper_N',FE.hNplus+E.DeltaN);
 A=[A;q]; %#ok<AGROW>
end
f=fopen(fullfile(here,'interval_scope.json'),'w');fprintf(f,'%s',jsonencode(A));fclose(f);
disp(A);fprintf('INTERVAL_SCOPE_DONE\n');
end
