function R=rev4_params()
%REV4_PARAMS Final numerical-validation settings.
%
% These settings parameterize numerical validation only.  They do not
% modify the controller, plant, actuator, or analytic theorem.

R.revision='Rev.4 Final Numerical Validation / Macro-Cascade Certificate';

% Local theorem-instance scenario.
R.certMode='certified_small';
R.certScale=0.125;
R.certReferenceTimes=0:5:120;
R.certRealizedTimes=30:5:120;
R.certHoldTime=10;

% Numerical differentiation.
R.relStep=5e-5;
R.stepCheck=[1e-4 5e-5 2.5e-5];

% Acceptance tolerances.
R.roundTripTol=5e-6;
R.phiTol=1e-8;
R.hurwitzTol=1e-8;
R.commonQTol=1e-10;
R.macroToWindTol=1e-8;
R.expectedFamilyCount=45;

% Macro architecture, primitive order (z,s,w,p,G,F,R,A).
R.macroPrimitiveIDs=[1 2 4 5 6 7 8];
R.windPrimitiveID=3;

% Results/export.
R.outputDir='rev4_outputs';
R.makePlots=true;
R.saveMat=true;

% Theory synchronization lock.
R.theoryFile='Port-Hamiltonian_UAV_IEEE_Theory_Sections_II-VIII_MacroCascade_Repaired_Rev2.docx';
R.theorySHA256='226cd1601a6617856dba6002a1a265d96e29552f3d5eac6cd3f50923340699a3';

% Interpretation locks.
R.claim.local='sampled numerical instantiation of the analytic local macro-cascade certificate';
R.claim.extended='nonlinear robustness/performance validation; not an Eq.-premise proof over the entire realized tube';
end
