function P=w2e_empty_pair()
%W2E_EMPTY_PAIR Exact paired-result schema.
P=struct('Index',0,'Proposed',w2e_empty_run(),'MatchedQP',w2e_empty_run(),'Pair',w2e_empty_pair_row());
end
