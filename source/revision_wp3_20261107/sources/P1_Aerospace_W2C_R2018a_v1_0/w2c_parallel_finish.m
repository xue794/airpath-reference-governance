function S=w2c_parallel_finish(pool,P,sourceRoot,packageRoot,token)
%W2C_PARALLEL_FINISH Restore worker paths; delete only runner-created pool.
S=struct('Restored',false,'DeletedCreatedPool',false,'Message','No pool.');if isempty(pool),return;end
try
 if ~isvalid(pool),return;end
 spmd,restored=w2c_worker_state('restore',sourceRoot,packageRoot,token);end %#ok<NASGU>
 S.Restored=true;S.Message='Worker paths restored.';
catch ME,S.Message=['Worker restore note: ' ME.identifier ': ' ME.message];end
if P.CreatedByRunner
 try,if isvalid(pool),delete(pool);end;S.DeletedCreatedPool=true;S.Message=[S.Message ' Runner-created pool closed.'];catch ME,S.Message=[S.Message ' Pool close note: ' ME.message];end
else,S.Message=[S.Message ' Pre-existing user pool retained.'];end
end
