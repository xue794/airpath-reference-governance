function [u,D] = w2f_qp_selector_core(u0,C,uB,O)
%W2F_QP_SELECTOR_CORE Matched-QP selector with prebuilt C and uB.
OO=O;OO.FeasibleStart=uB(:);
[u,D]=w2b_qp3(u0(:),-C.G,C.c,O.Scale,OO);
if ~D.Success
    error('W2F:QP','No verified KKT solution: %s',D.Status);
end
end
