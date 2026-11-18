function [P,pool]=w2d_parallel_prepare(useParallel,requested,numJobs,sourceRoot,packageRoot,outdir)
%W2D_PARALLEL_PREPARE Local process pool capped to independent trajectory count.
pool=[];effective=min(requested,numJobs);P=struct('Requested',useParallel,'UserRequestedWorkers',requested,'EffectiveRequestedWorkers',effective,'Enabled',false,'Workers',0,'PoolSize',0,'CreatedByRunner',false,'ReusedPool',false,'Mode','SERIAL','Message','','WorkerInfo',{{}});
if ~useParallel||effective<2,P.Message='Serial deterministic trajectory execution.';return;end
try,available=exist('parpool','file')==2 && license('test','Distrib_Computing_Toolbox');catch,available=false;end
if ~available,P.Message='Parallel Computing Toolbox unavailable; using serial trajectories.';fprintf('%s\n',P.Message);return;end
try
 pool=gcp('nocreate');if ~isempty(pool),P.ReusedPool=true;else
  c=parcluster('local');if c.NumWorkers<effective,c.NumWorkers=effective;end
  fprintf('Starting local PROCESS pool with %d workers (W2-D has %d independent trajectories; user requested %d).\n',effective,numJobs,requested);
  pool=parpool(c,effective);P.CreatedByRunner=true;
 end
 if isempty(pool),error('W2D:NoPool','No pool available.');end
 token=outdir;spmd,workerState=w2d_worker_state('initialize',sourceRoot,packageRoot,token);end
 W=cell(pool.NumWorkers,1);for k=1:pool.NumWorkers,W{k}=workerState{k};end;P.WorkerInfo=W;P.PoolSize=pool.NumWorkers;P.Enabled=true;P.Workers=min(effective,pool.NumWorkers);P.Mode='LOCAL_PROCESS_PARFOR';P.Message='Independent trajectories parallel; reporting/timing client-serial.';
 try,writetable(struct2table([W{:}]),fullfile(outdir,'W2D_worker_identity.csv'));catch,end
 fprintf('Parallel enabled: %d workers used for %d independent trajectories.\n',P.Workers,numJobs);
catch ME
 P.Message=['Parallel setup unavailable: ' ME.identifier ': ' ME.message ' | serial fallback.'];fprintf(2,'%s\n',P.Message);
 try,if ~isempty(pool)&&P.CreatedByRunner,delete(pool);end;catch,end;pool=[];P.Enabled=false;P.Workers=0;P.Mode='SERIAL_SETUP_FALLBACK';
end
end
