function C=w2e_compare_samples(A,B)
%W2E_COMPARE_SAMPLES Exact field-by-field comparison of complete uncertainty records.
C=struct('ExactMatch',false,'N',numel(A),'MismatchCount',0,'FirstMismatchIndex',NaN,'FirstMismatchField','');
if numel(A)~=numel(B),C.MismatchCount=abs(numel(A)-numel(B))+1;return;end
fa=fieldnames(A);fb=fieldnames(B);if ~isequal(fa,fb),C.MismatchCount=1;C.FirstMismatchField='FIELD_SET_OR_ORDER';return;end
for k=1:numel(A)
 for j=1:numel(fa)
  f=fa{j};if ~isequaln(A(k).(f),B(k).(f)),C.MismatchCount=C.MismatchCount+1;if isnan(C.FirstMismatchIndex),C.FirstMismatchIndex=k;C.FirstMismatchField=f;end;end
 end
end
C.ExactMatch=C.MismatchCount==0;
end
