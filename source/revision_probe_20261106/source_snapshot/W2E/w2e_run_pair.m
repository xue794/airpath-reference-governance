function P=w2e_run_pair(Samp,Rmc,Cbase,O)
%W2E_RUN_PAIR Same frozen uncertainty record is used for both selectors.
% The two method runs are independently exception-contained: one failure never suppresses the counterpart.
A=w2e_run_method(Samp,Rmc,Cbase,'Proposed',O);
B=w2e_run_method(Samp,Rmc,Cbase,'MatchedQP',O);
P=w2e_empty_pair();P.Index=Samp.index;P.Proposed=A;P.MatchedQP=B;P.Pair=w2e_pair_row(A,B);
end
