function P=w2c_pairwise_row(SP,SQ,EP,EQ)
%W2C_PAIRWISE_ROW Same-mode summary. Own-trajectory metrics are not pointwise optimality comparisons.
assert(strcmp(SP.mode,SQ.mode),'W2C:Pair','Modes differ.');
P=struct();P.Mode=SP.mode;P.InitialStateInfDiff=max(abs(SP.x(1,:)-SQ.x(1,:)));
P.QP_minus_Proposed_RMS_m=EQ.PosRMS_Full_m-EP.PosRMS_Full_m;
P.QP_minus_Proposed_Peak_m=EQ.PosPeak_Full_m-EP.PosPeak_Full_m;
P.QP_minus_Proposed_LatePeak_m=EQ.PosPeak_Late_m-EP.PosPeak_Late_m;
P.QP_minus_Proposed_Intervention_pp=EQ.GovInterventionPct-EP.GovInterventionPct;
if EP.PosRMS_Full_m>0,P.RMSChangePct=100*(EQ.PosRMS_Full_m/EP.PosRMS_Full_m-1);else,P.RMSChangePct=NaN;end
if EP.PosPeak_Full_m>0,P.PeakChangePct=100*(EQ.PosPeak_Full_m/EP.PosPeak_Full_m-1);else,P.PeakChangePct=NaN;end
P.MaxStateGridInfDiff=max(max(abs(SP.x-SQ.x)));
physicalRows=11:17;thrustRows=[16 17];
P.QPNonbackupPhysicalSupportSamples=sum(any(SQ.Y.qpSupport(:,physicalRows),2));
P.QPThrustSupportSamples=sum(any(SQ.Y.qpSupport(:,thrustRows),2));
P.QPTminSupportSamples=sum(SQ.Y.qpSupport(:,16));P.QPTmaxSupportSamples=sum(SQ.Y.qpSupport(:,17));
P.QPNonbackupPhysicalActiveSetSamples=sum(any(SQ.Y.qpActive(:,physicalRows),2));
P.QPThrustActiveSetSamples=sum(any(SQ.Y.qpActive(:,thrustRows),2));
P.QPReferenceGatePass=EQ.ReferenceGatePass;P.QPKKTGatePass=EQ.QPKKTGatePass;P.ProposedReferenceGatePass=EP.ReferenceGatePass;
end
