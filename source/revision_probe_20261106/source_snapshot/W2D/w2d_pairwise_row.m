function P=w2d_pairwise_row(SP,SQ,EP,EQ)
%W2D_PAIRWISE_ROW Same targeted mode, own-trajectory comparison and evidence counts.
assert(strcmp(SP.mode,SQ.mode),'W2D:Pair','Modes differ.');
P=struct();P.Mode=SP.mode;P.ExpectedRow=SP.Case.ExpectedRow;P.ExpectedRowIndex=SP.Case.ExpectedRowIndex;P.InitialStateInfDiff=max(abs(SP.x(1,:)-SQ.x(1,:)));
P.QP_minus_Proposed_RMS_m=EQ.PosRMS_Full_m-EP.PosRMS_Full_m;P.QP_minus_Proposed_Peak_m=EQ.PosPeak_Full_m-EP.PosPeak_Full_m;P.QP_minus_Proposed_LatePeak_m=EQ.PosPeak_Late_m-EP.PosPeak_Late_m;P.QP_minus_Proposed_Intervention_pp=EQ.GovInterventionPct-EP.GovInterventionPct;
if EP.PosRMS_Full_m>0,P.RMSChangePct=100*(EQ.PosRMS_Full_m/EP.PosRMS_Full_m-1);else,P.RMSChangePct=NaN;end
if EP.PosPeak_Full_m>0,P.PeakChangePct=100*(EQ.PosPeak_Full_m/EP.PosPeak_Full_m-1);else,P.PeakChangePct=NaN;end
P.MaxStateGridInfDiff=max(max(abs(SP.x-SQ.x)));idx=SP.Case.ExpectedRowIndex;
P.ProposedExpectedLimiterSamples=sum(SP.Y.scalarLimiter(:,idx));P.QPExpectedSupportSamples=sum(SQ.Y.qpSupport(:,idx));P.QPExpectedActiveSetSamples=sum(SQ.Y.qpActive(:,idx));P.QPExpectedBindingSamples=sum(SQ.Y.binding(:,idx));
P.QPNonbackupPhysicalSupportSamples=sum(any(SQ.Y.qpSupport(:,11:17),2));P.QPBackupSupportSamples=sum(SQ.Y.qpSupport(:,18));
ix=find(SQ.Y.qpSupport(:,idx),1,'first');if isempty(ix),P.FirstQPExpectedSupportTime_s=NaN;else,P.FirstQPExpectedSupportTime_s=SQ.t(ix);end
P.QPReferenceGatePass=EQ.ReferenceGatePass;P.QPKKTGatePass=EQ.QPKKTGatePass;P.ProposedReferenceGatePass=EP.ReferenceGatePass;
end
