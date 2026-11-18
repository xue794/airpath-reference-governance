function D=viii_instance_data_rev2(C,P8)
%VIII_INSTANCE_DATA_REV2 Reuse the Rev.1.1.2 trajectory cache when available.

if nargin<1 || isempty(C), C=cl6_params(); end
if nargin<2 || isempty(P8), P8=viii_params(); end

candidates={ ...
    fullfile(pwd,'viii_preflight_cache_rev112.mat'), ...
    fullfile(fileparts(pwd), ...
      'T2_GTM_GNA_SectionVIII_CompositeCertificate_Rev1_1_2', ...
      'viii_preflight_cache_rev112.mat')};

for k=1:numel(candidates)
    if exist(candidates{k},'file')==2
        S=load(candidates{k},'D');
        if isfield(S,'D') && isfield(S.D,'varying') && isfield(S.D,'noisy')
            fprintf('Loading accepted varying/noisy trajectories from:\n   %s\n',candidates{k});
            D=S.D;
            return;
        end
    end
end

fprintf('No compatible Rev.1.1.2 cache found; generating the two accepted cases once.\n');
D.varying=simulate_complete_6dof_case('varying',C);
D.noisy=simulate_complete_6dof_case('noisy',C);

try
    save(fullfile(pwd,'viii_preflight_cache_rev112.mat'),'D','-v7');
catch ME
    warning('viii_instance_data_rev2:cache','Could not save cache: %s',ME.message);
end
end
