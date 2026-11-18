function E=w2e_run_method(Samp,Rmc,Cbase,method,O)
%W2E_RUN_METHOD Exception-contained execution of one method on one frozen sample.
E=w2e_empty_run();E.Index=Samp.index;E.Method=method;E.rhoB=O.RhoB;
try
 C=rev41_apply_sample(Cbase,Samp,Rmc);assert(abs(C.R.rhoB-O.RhoB)<=1e-14,'W2E:RhoRuntime','rhoB changed during sample application.');
 S=w2e_simulate_trial(method,C,O);E=w2e_summarize_run(S,Samp,Rmc,O);
catch ME
 E.Status='FAILED';E.Completed=false;E.ErrorId=ME.identifier;E.ErrorMessage=ME.message;[E.FailureClass,E.DomainExit,E.NumericalFailure]=w2e_classify_error(ME.identifier);
end
end
