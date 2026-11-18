function P=viii_collect_witness_points(Data,C,P8)
%VIII_COLLECT_WITNESS_POINTS Feasible centers for the Rev.1.1 macro witness.
%
% Data.varying and Data.noisy are accepted complete-6DOF trajectories.

P=struct('t',{},'X',{},'mode',{},'kind',{});
n=0;
B=viii_blocks();
Z=zeros(B.n,1);

% Straight/hold zero centers are genuinely feasible.
for kt=1:numel(P8.sample.hold)
    n=n+1;
    P(n).t=P8.sample.hold(kt);
    P(n).X=Z;
    P(n).mode='hold';
    P(n).kind='zero-hold';
end

cases={Data.varying,Data.noisy};
ts={P8.sample.actualVarying,P8.sample.actualNoisy};
names={'varying','noisy'};

for ic=1:2
    S=cases{ic};
    for kt=1:numel(ts{ic})
        [~,k]=min(abs(S.t-ts{ic}(kt)));
        [X,~]=viii_reduced_from_full(S.t(k),S.x(k,:).',names{ic},C,P8);

        n=n+1;
        P(n).t=S.t(k);
        P(n).X=X;
        P(n).mode=names{ic};
        P(n).kind='realized';
    end
end
end
