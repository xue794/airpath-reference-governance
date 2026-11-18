function O=w2d_options()
%W2D_OPTIONS Frozen, predeclared targeted constraint-activation protocol.
O.Version='W2D_R2018a_v1.0';
O.RhoB=0.10;
O.Methods={'Proposed','MatchedQP'};
O.Modes={'target_tmin','target_tmax','target_nplus'};
O.UseParallel=true;
O.RequestedWorkers=6;
O.SaveFullTrajectories=true;
% The target set is inherited from the W2-B finite witness scan.  It is not
% tuned after observing W2-C closed-loop results.
O.RampStart_s=5.0;
O.RampDuration_s=0.50;      % C2 quintic rise/fall.
O.PulseHold_s=0.75;         % Exact target velocity is held briefly.
O.CaseDuration_s=20.0;
O.LateStart_s=12.0;
O.TargetTmin=[30.0;0;0];
O.TargetTmax=[46.0;0;0];
O.TargetNplus=[38.1;8*pi/180;0];
O.InterventionTol=1e-10;
O.LimiterTol=1e-10;
O.BindingRelTol=1e-8;
O.NormalizedPhiTol=1e-8;
O.BackupTol=1e-8;
O.KKTGateTol=5e-8;
O.InitialBackupTol=1e-12;
O.QP=w2b_options();
assert(abs(O.QP.RhoB-O.RhoB)<=1e-14,'W2D:Rho','W2-B QP protocol rho mismatch.');
assert(norm(O.QP.Scale-[2;8*pi/180;12*pi/180],inf)<1e-14,'W2D:Scale','Locked QP scale changed.');
end
