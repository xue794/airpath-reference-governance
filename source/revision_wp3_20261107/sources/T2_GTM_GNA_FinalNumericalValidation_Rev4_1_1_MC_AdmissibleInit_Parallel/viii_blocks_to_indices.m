function idx=viii_blocks_to_indices(blockIds)
%VIII_BLOCKS_TO_INDICES Convert atomic block IDs to reduced-coordinate indices.
B=viii_blocks();
idx=[];
for k=1:numel(blockIds)
    idx=[idx B.idx{blockIds(k)}]; %#ok<AGROW>
end
idx=idx(:).';
end
