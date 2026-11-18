function run_representative_direct(outputRoot)
% Reintegrate the archived C1 Direct case in an isolated output tree.
% MATLAB R2018a. No model, controller, parameter or integration changes.
releaseRoot=fileparts(fileparts(mfilename('fullpath')));
if nargin<1
    outputRoot=fullfile(releaseRoot,'generated','matlab_C1_005');
end
assert(ischar(outputRoot),'Supply a character-vector output path.');
assert(~exist(outputRoot,'dir') && ~exist(outputRoot,'file'), ...
    'Output already exists. Inspect it and choose a new directory; no run is overwritten.');
assert(exist(fullfile(releaseRoot,'source','decel_center_20261108','run_direct.m'),'file')==2);
[ok,msg]=mkdir(outputRoot);assert(ok,msg);
% Resolve a caller-relative output before changing to the copied study folder.
% fileattrib is available in the retained MATLAB R2018a environment.
[ok,outputAttributes]=fileattrib(outputRoot);assert(ok,'Cannot resolve output directory.');
assert(isscalar(outputAttributes) && outputAttributes.directory,'Expected one output directory.');
outputRoot=outputAttributes.Name;
[ok,msg]=copyfile(fullfile(releaseRoot,'source'),fullfile(outputRoot,'source'));assert(ok,msg);
oldPath=path;oldFolder=pwd;
cleanup=onCleanup(@() restoreEnvironment(oldPath,oldFolder)); %#ok<NASGU>
directRoot=fullfile(outputRoot,'source','decel_center_20261108');
cd(directRoot);addpath(directRoot);
started=tic;
run_direct('C1',.005);
record=struct('case','C1','arm','Direct','max_step_s',.005, ...
    'matlab_version',version,'total_wall_seconds',toc(started), ...
    'source_policy','Unmodified archived source copied to isolated output tree');
fid=fopen(fullfile(outputRoot,'release_run.json'),'w');assert(fid~=-1);
fprintf(fid,'%s',jsonencode(record));fclose(fid);
fprintf('REPRESENTATIVE_RUN_DONE %s\n',outputRoot);
end

function restoreEnvironment(savedPath,savedFolder)
% Anonymous capture preserves values after the main workspace exits in R2018a.
path(savedPath);cd(savedFolder);
end
