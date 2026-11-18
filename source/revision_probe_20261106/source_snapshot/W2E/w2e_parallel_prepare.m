function [P,pool]=w2e_parallel_prepare(useParallel,requested,numPending,sourceRoot,packageRoot,outdir)
%W2E_PARALLEL_PREPARE Pair-level process parallelism; each worker executes both methods for one sample.
pool=[];effective=min(requested,max(1,numPending));P=struct('Requested',useParallel,'UserRequestedWorkers',requested,'EffectiveRequestedWorkers',effective,'Enabled',false,'Workers',0,'PoolSize',0,'CreatedByRunner',false,'ReusedPool',false,'Mode','SERIAL','Message','','WorkerInfo',{{}});
if ~useParallel||effective<2,P.Message='Serial paired Monte-Carlo execution.';return;end
try,available=exist('parpool','file')==2&&license('test','Distrib_Computing_Toolbox');catch,available=false;end
if ~available,P.Message='Parallel Computing Toolbox unavailable; using serial execution.';fprintf('%s\n',P.Message);return;end
try
 pool=gcp('nocreate');if ~isempty(pool),P.ReusedPool=true;if pool.NumWorkers<effective,fprintf(2,'Pre-existing pool has %d workers; requested %d. It is retained and used without deletion.\n',pool.NumWorkers,effective);end;else,c=parcluster('local');if c.NumWorkers<effective,c.NumWorkers=effective;end;fprintf('Starting local PROCESS pool with %d workers for pair-level Monte Carlo...\n',effective);pool=parpool(c,effective);P.CreatedByRunner=true;end
 if isempty(pool),error('W2E:NoPool','No pool available.');end
 token=outdir;spmd,workerState=w2e_worker_state('initialize',sourceRoot,packageRoot,token);end
 W=cell(pool.NumWorkers,1);for k=1:pool.NumWorkers,W{k}=workerState{k};end;P.WorkerInfo=W;P.PoolSize=pool.NumWorkers;P.Enabled=true;P.Workers=min(numPending,pool.NumWorkers);P.Mode='LOCAL_PROCESS_PAIR_PARFOR';P.Message='Each PARFOR iteration runs one complete Proposed/Matched-QP pair on the same saved sample.';try,writetable(struct2table([W{:}]),fullfile(outdir,'W2E_worker_identity.csv'));catch,end;fprintf('Parallel enabled: pool=%d workers, pending pairs=%d.\n',pool.NumWorkers,numPending);
catch ME,P.Message=['Parallel setup unavailable: ' ME.identifier ': ' ME.message ' | serial fallback.'];fprintf(2,'%s\n',P.Message);try,if ~isempty(pool)&&P.CreatedByRunner,delete(pool);end;catch,end;pool=[];P.Enabled=false;P.Workers=0;P.Mode='SERIAL_SETUP_FALLBACK';
end
end
