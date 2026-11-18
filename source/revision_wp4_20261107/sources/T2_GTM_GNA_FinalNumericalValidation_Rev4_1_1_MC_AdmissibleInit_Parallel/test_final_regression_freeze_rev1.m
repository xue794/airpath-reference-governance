function Results=test_final_regression_freeze_rev1()
%TEST_FINAL_REGRESSION_FREEZE_REV1 Final one-command acceptance regression.
%
% The complete Rev.1 test recursively runs the verified Section-VII chain:
% VII -> VI-D -> VI-C -> VI-B -> VI-A -> V -> IV -> plant.
% Thus the synchronized parameter set is exercised through every layer once.

Sync=test_final_parameter_sync_rev1();
CL=test_complete_6dof_rev1();

fprintf('\nFINAL REGRESSION FREEZE REV.1 ACCEPTANCE\n');
fprintf('========================================\n');
fprintf('All synchronized module and complete-6DOF assertions passed.\n');
fprintf('Canonical K_omega,f = %.1f I.\n',CL.C.F.KomegaF(1,1));
fprintf('No integration-only controller override remains.\n');

Results=struct('Sync',Sync,'CL',CL);
end
