function P=viii_macro_partitions()
%VIII_MACRO_PARTITIONS Hierarchical two-macro-subsystem candidate partitions.
%
% Atomic block order in the project is (z,s,w,p,G,F,R,A).
% Candidate partitions preserve the architecture and progressively absorb
% strongly coupled inner blocks into the outer macro subsystem.

P=struct('name',{},'g1Blocks',{},'g2Blocks',{});

P(1).name='Outer(z,s,w,G) | Inner(p,F,R,A)';
P(1).g1Blocks=[1 2 3 5];
P(1).g2Blocks=[4 6 7 8];

P(2).name='Outer+Path(z,s,w,G,p) | Attitude+Actuator(F,R,A)';
P(2).g1Blocks=[1 2 3 5 4];
P(2).g2Blocks=[6 7 8];

P(3).name='Outer+Path+Filter | Rotation+Actuator';
P(3).g1Blocks=[1 2 3 5 4 6];
P(3).g2Blocks=[7 8];

P(4).name='All controlled states except actuator | Actuator';
P(4).g1Blocks=[1 2 3 5 4 6 7];
P(4).g2Blocks=8;
end
