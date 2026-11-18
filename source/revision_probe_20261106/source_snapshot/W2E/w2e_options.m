function O=w2e_options()
%W2E_OPTIONS Frozen paired Monte-Carlo protocol for Paper-1 Aerospace revision.
% The protocol is fixed before seeing any W2-E outcomes.
root=fileparts(mfilename('fullpath'));
O.Version='W2E_R2018a_v1.0';
O.ProtocolID='P1_Aerospace_W2E_PairedMC500_rho010_v1_0';
O.RhoB=0.10;
O.N=500;
O.Methods={'Proposed','MatchedQP'};
O.Mode='monte_carlo';
O.UseParallel=true;
O.RequestedWorkers=56;
O.BatchPairsPerWorker=2; % 2 worker-waves per checkpoint, as in frozen Rev.4.1.1 parallel campaign.
O.Resume=true;
O.MakePlots=false;       % publication figures are frozen only after audit.
O.SampleSource=fullfile(root,'data','Rev41Results_FROZEN_SAMPLE_SOURCE.mat');
O.ExpectedOldSeed=20260819;
O.ExpectedOldN=500;
% Selector/reference diagnostics. These do NOT soften any constraint.
O.InterventionTol=1e-10;
O.LimiterTol=1e-10;
O.BindingRelTol=1e-8;
O.NormalizedPhiTol=1e-8;
O.BackupTol=1e-8;
O.KKTGateTol=5e-8;
% Primary paired continuous-performance population: BOTH methods must finish
% and pass the frozen CorePass gate. No failure penalty is imputed.
O.PairedPopulation='BOTH_CORE_PASS';
O.QP=w2b_options();
assert(abs(O.QP.RhoB-O.RhoB)<=1e-14,'W2E:Rho','Frozen QP rho does not match W2-E.');
assert(norm(O.QP.Scale-[2;8*pi/180;12*pi/180],inf)<1e-14,'W2E:Scale','Frozen QP metric changed.');
end
