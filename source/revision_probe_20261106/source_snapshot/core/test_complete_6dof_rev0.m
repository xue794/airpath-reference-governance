function Results=test_complete_6dof_rev0()
%TEST_COMPLETE_6DOF_REV0 Compatibility wrapper.
% The canonical controller parameter set has been synchronized and frozen
% in Final Regression Freeze Rev.1.
warning('test_complete_6dof_rev0:Redirect', ...
    'Rev.0 entry point redirects to final synchronized Rev.1.');
Results=test_complete_6dof_rev1();
end
