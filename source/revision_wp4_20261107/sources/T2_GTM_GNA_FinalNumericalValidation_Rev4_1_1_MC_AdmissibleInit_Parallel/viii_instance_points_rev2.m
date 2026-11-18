function P=viii_instance_points_rev2(Data,C,P8,P2)
%VIII_INSTANCE_POINTS_REV2 Feasible centers for numerical theorem instantiation.

P=struct('t',{},'X',{},'mode',{},'kind',{});
n=0;
B=viii_blocks();
Z=zeros(B.n,1);

% Straight hold: zero reduced state is feasible.
for kt=1:numel(P2.sample.hold)
    n=n+1;
    P(n).t=P2.sample.hold(kt);
    P(n).X=Z;
    P(n).mode='hold';
    P(n).kind='zero-hold';
end

cases={Data.varying,Data.noisy};
ts={P2.sample.varying,P2.sample.noisy};
names={'varying','noisy'};

for ic=1:2
    S=cases{ic};
    for kt=1:numel(ts{ic})
        [~,k]=min(abs(S.t-ts{ic}(kt)));
        X=viii_reduced_from_full(S.t(k),S.x(k,:).',names{ic},C,P8);
        n=n+1;
        P(n).t=S.t(k);
        P(n).X=X;
        P(n).mode=names{ic};
        P(n).kind='realized';
    end
end
end
