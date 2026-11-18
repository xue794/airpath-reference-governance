function Results=run_complete_6dof_rev0(makeFigures)
%RUN_COMPLETE_6DOF_REV0 Compatibility wrapper for final synchronized Rev.1.
if nargin<1, makeFigures=false; end
warning('run_complete_6dof_rev0:Redirect', ...
    'Rev.0 entry point redirects to final synchronized Rev.1.');
Results=run_complete_6dof_rev1(makeFigures);
end
