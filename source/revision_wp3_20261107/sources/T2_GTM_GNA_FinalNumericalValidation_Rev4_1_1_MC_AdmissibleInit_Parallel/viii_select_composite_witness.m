function W=viii_select_composite_witness(Data,C,P8)
%VIII_SELECT_COMPOSITE_WITNESS Prefer macro small gain, otherwise direct composite.
%
% The formal manuscript theorem remains small-gain/composite and condition-based.
% Numerically, a direct full-state common quadratic witness is accepted as a
% stronger sampled local composite certificate if the macro split is not suitable.

M=viii_macro_small_gain_witness(Data,C,P8);

if M.valid
    W=M;
    W.mode='macro';
    W.selectionMessage='two-macro small-gain witness';
    return;
end

D=viii_direct_composite_witness(Data,C,P8,M.J);
D.macroDiagnostic=M;

if D.valid
    W=D;
    W.selectionMessage='direct 34-state sampled composite quadratic witness';
else
    W=D;
    W.selectionMessage='no sampled composite witness closed';
end
end
