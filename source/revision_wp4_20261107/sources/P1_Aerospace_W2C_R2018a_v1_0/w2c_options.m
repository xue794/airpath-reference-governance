function O=w2c_options()
%W2C_OPTIONS Frozen deterministic comparison protocol.
O.Version='W2C_R2018a_v1.0';
O.RhoB=0.10;
O.Methods={'Proposed','MatchedQP'};
O.Modes={'hold','varying','noisy'};
O.UseParallel=true;
O.RequestedWorkers=6;       % Runner caps to number of independent trajectories.
O.SaveFullTrajectories=true;
O.InterventionTol=1e-10;
O.LimiterTol=1e-10;
O.BindingRelTol=1e-8;
O.NormalizedPhiTol=1e-8;
O.BackupTol=1e-8;
O.HoldStateMatchTol=1e-7;
O.KKTGateTol=5e-8;
O.TimingEnabled=true;
O.TimingSamplesPerStressCase=40;
O.TimingWarmup=3;
O.TimingRepeats=10;
O.QP=w2b_options();
assert(abs(O.QP.RhoB-O.RhoB)<=1e-14,'W2C:Rho','W2-B QP protocol rho mismatch.');
end
