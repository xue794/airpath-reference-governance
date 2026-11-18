function J=w2c_run_job(mode,method,O)
%W2C_RUN_JOB Exception-contained independent trajectory job.
J=struct('Success',false,'Mode',mode,'Method',method,'S',[],'Entry',[], ...
 'ErrorId','','ErrorMessage','','ErrorStack','');
try
 C=cl6_params();assert(abs(C.R.rhoB-.20)<1e-12,'W2C:Baseline','Expected untouched source default rhoB=0.20.');
 C.R.rhoB=O.RhoB;
 S=w2c_simulate_case(mode,method,C,O);
 E=w2c_summarize_case(S,O);
 J.Success=true;J.S=S;J.Entry=E;
catch ME
 J.ErrorId=ME.identifier;J.ErrorMessage=ME.message;
 st='';for k=1:numel(ME.stack),st=[st sprintf('#%d %s line %d\n',k,ME.stack(k).name,ME.stack(k).line)];end %#ok<AGROW>
 J.ErrorStack=st;
 E=struct('Method',method,'Mode',mode,'Status','FAILED','ErrorId',ME.identifier,'ErrorMessage',ME.message,'rhoB',O.RhoB);
 J.Entry=E;
end
end
