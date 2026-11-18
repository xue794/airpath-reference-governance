function D=viii_preflight_data(C,P8)
%VIII_PREFLIGHT_DATA Generate or load the two complete trajectories needed by Section VIII.
%
% Rev.1.1.2 caches the expensive varying/noisy trajectories in the package
% directory.  If a later certificate-only change is needed, the ODE cases do
% not have to be recomputed.

if nargin<2 || isempty(P8)
    P8=viii_params();
end

cacheFile=fullfile(pwd,P8.cache.file);

if P8.cache.enable && exist(cacheFile,'file')==2
    S=load(cacheFile,'D','signature');
    sigNow=local_signature(C,P8);
    if isfield(S,'signature') && isequal(S.signature,sigNow) && isfield(S,'D')
        fprintf('\nLoading Section-VIII preflight trajectories from cache:\n');
        fprintf('   %s\n',cacheFile);
        D=S.D;
        return;
    end
end

fprintf('\nGenerating Section-VIII preflight complete-loop trajectories...\n');
fprintf('   varying wind case...\n');
D.varying=simulate_complete_6dof_case('varying',C);
fprintf('   bounded-measurement-error case...\n');
D.noisy=simulate_complete_6dof_case('noisy',C);

if P8.cache.enable
    signature=local_signature(C,P8); %#ok<NASGU>
    try
        save(cacheFile,'D','signature','-v7');
        fprintf('   cached trajectories -> %s\n',cacheFile);
    catch ME
        warning('viii_preflight_data:CacheWrite', ...
            'Could not write preflight cache: %s',ME.message);
    end
end
end

function sig=local_signature(C,P8)
sig=[C.F.kRf;diag(C.F.KomegaF);C.Dp.kR;diag(C.Dp.KomegaC); ...
     diag(C.G.Lambda);diag(C.G.Ks);C.G.Kc(:);C.E.Lw(1,1); ...
     P8.Vstar;P8.flowDt;P8.jacStep];
end
