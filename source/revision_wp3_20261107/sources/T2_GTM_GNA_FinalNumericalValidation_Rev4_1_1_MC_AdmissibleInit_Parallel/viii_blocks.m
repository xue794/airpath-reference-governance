function B=viii_blocks()
%VIII_BLOCKS Reduced-state block indices, order (z,s,w,p,G,F,R,A).
B.idx={1:3,4:6,7:9,10:12,13:18,19:24,25:30,31:34};
B.dim=[3 3 3 3 6 6 6 4];
B.n=34;
B.labels={'z','s','w','p','G','F','R','A'};
end
