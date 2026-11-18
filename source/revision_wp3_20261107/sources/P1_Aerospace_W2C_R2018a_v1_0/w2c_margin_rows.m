function [P,M]=w2c_margin_rows(S)
%W2C_MARGIN_ROWS Per-affine-row minima and state-level feasibility margins.
P=cell(18,1);for j=1:18,P{j}=struct('Method',S.Method,'Mode',S.mode,'RowIndex',j,'Row',S.Y.constraintNames{j},'MinAppliedPhi',min(S.Y.phi(:,j)),'MinNominalPhi',min(S.Y.phi0(:,j)));end
M=cell(numel(S.Y.stateMarginNames),1);for j=1:numel(M),M{j}=struct('Method',S.Method,'Mode',S.mode,'Margin',S.Y.stateMarginNames{j},'MinValue',min(S.Y.stateMargins(:,j)));end
end
