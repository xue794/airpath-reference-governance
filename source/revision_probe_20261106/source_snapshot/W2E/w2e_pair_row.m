function P=w2e_pair_row(A,B)
%W2E_PAIR_ROW Paired outcomes. Continuous deltas exist ONLY on BothCorePass population.
% Delta convention is always MatchedQP minus Proposed.
P=w2e_empty_pair_row();P.Index=A.Index;assert(A.Index==B.Index,'W2E:PairIndex','Pair indices differ.');P.ProposedStatus=A.Status;P.QPStatus=B.Status;P.ProposedFailureClass=A.FailureClass;P.QPFailureClass=B.FailureClass;
P.BothCompleted=A.Completed&&B.Completed;P.BothHardFeasible=A.HardFeasible&&B.HardFeasible;P.BothCorePass=A.CorePass&&B.CorePass;P.BothCertifiedDomain=A.CertifiedDomain&&B.CertifiedDomain;P.ProposedOnlyCorePass=A.CorePass&&~B.CorePass;P.QPOnlyCorePass=B.CorePass&&~A.CorePass;P.NeitherCorePass=~A.CorePass&&~B.CorePass;P.ExecutionDiscordant=xor(A.Completed,B.Completed);P.PairPerformanceEligible=P.BothCorePass;
if P.BothCorePass
 P.OutcomeCode='BOTH_CORE_PASS';
 P.dPosRMSFull_m=B.PosRMSFull_m-A.PosRMSFull_m;P.dPosPeakFull_m=B.PosPeakFull_m-A.PosPeakFull_m;P.dPosRMSLate_m=B.PosRMSLate_m-A.PosRMSLate_m;P.dPosPeakLate_m=B.PosPeakLate_m-A.PosPeakLate_m;P.dIntervention_pp=B.GovInterventionPct-A.GovInterventionPct;P.dMeanCorrectionNorm=B.MeanCorrectionNorm-A.MeanCorrectionNorm;P.dIntegratedObjective_s=B.IntegratedObjective_s-A.IntegratedObjective_s;P.dMaxAttitudeTube_deg=B.MaxAttitudeTube_deg-A.MaxAttitudeTube_deg;P.dMaxMomentResidual_Nm=B.MaxMomentResidual_Nm-A.MaxMomentResidual_Nm;P.dMaxThrustRate_Nps=B.MaxThrustRate_Nps-A.MaxThrustRate_Nps;
elseif P.ProposedOnlyCorePass,P.OutcomeCode='PROPOSED_ONLY_CORE_PASS';
elseif P.QPOnlyCorePass,P.OutcomeCode='QP_ONLY_CORE_PASS';
elseif P.BothCompleted,P.OutcomeCode='BOTH_COMPLETED_CORE_GATE_FAIL';
elseif P.ExecutionDiscordant,P.OutcomeCode='DISCORDANT_EXECUTION';
else,P.OutcomeCode='BOTH_NONCOMPLETE';end
end
