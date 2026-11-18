function S=w2e_worker_state(action,sourceRoot,packageRoot,token)
%W2E_WORKER_STATE Temporary worker paths and exact source/package verification.
persistent previousPath previousPwd activeToken
S=struct('Pass',true,'Action',action,'TaskId',NaN,'PackageRoot',packageRoot,'SourceRoot',sourceRoot,'QPPath','','NumComputationThreads',NaN);
if strcmp(action,'restore')
 if ~isempty(activeToken) && strcmp(activeToken,token)
  if ~isempty(previousPwd)&&exist(previousPwd,'dir')==7,cd(previousPwd);end;if ~isempty(previousPath),path(previousPath);end
  previousPath=[];previousPwd=[];activeToken=[];
 end;return;
end
assert(strcmp(action,'initialize'),'W2E:WorkerAction','Unknown worker action.');
if isempty(activeToken),previousPath=path;previousPwd=pwd;activeToken=token;elseif ~strcmp(activeToken,token),error('W2E:WorkerState','Worker has an unfinished W2-E session.');end
addpath(sourceRoot,'-begin');addpath(packageRoot,'-begin');cd(packageRoot);rehash path;
E=load(fullfile(packageRoot,'data','w2e_expected_source.mat'));
for k=1:numel(E.Names)
 f=fullfile(sourceRoot,E.Names{k});[~,fn]=fileparts(f);a=strrep(which(fn),char(92),'/');b=strrep(f,char(92),'/');if ispc,same=strcmpi(a,b);else,same=strcmp(a,b);end
 assert(same,'W2E:WorkerPath','Worker source shadow: %s -> %s',fn,a);
 actual=fileread(f);actual=strrep(actual,sprintf('\r\n'),sprintf('\n'));actual=strrep(actual,sprintf('\r'),sprintf('\n'));
 expected=E.Text{k};expected=strrep(expected,sprintf('\r\n'),sprintf('\n'));expected=strrep(expected,sprintf('\r'),sprintf('\n'));assert(strcmp(actual,expected),'W2E:WorkerSource','Worker source mismatch: %s',fn);
end
S.QPPath=which('w2b_qp3');a=strrep(S.QPPath,char(92),'/');b=strrep(fullfile(packageRoot,'w2b_qp3.m'),char(92),'/');if ispc,same=strcmpi(a,b);else,same=strcmp(a,b);end;assert(same,'W2E:WorkerPackage','Wrong W2-B solver on worker.');
try,t=getCurrentTask();if ~isempty(t),S.TaskId=t.ID;end;catch,end;try,S.NumComputationThreads=maxNumCompThreads;catch,end
end
