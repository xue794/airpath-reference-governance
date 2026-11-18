function O=viii_common_group_storage(Acells,P8)
%VIII_COMMON_GROUP_STORAGE Find a sampled common quadratic storage for one macro group.
%
% Base MATLAB only.  Candidate Lyapunov matrices are generated from every
% Hurwitz sampled self matrix and from the mean self matrix.  A deterministic
% convex-blending search then minimizes the worst generalized derivative rate.

nC=numel(Acells);
n=size(Acells{1},1);

alphaSpec=-inf(nC,1);
for k=1:nC
    alphaSpec(k)=max(real(eig(Acells{k})));
end

Pcand={};

Abar=zeros(n);
for k=1:nC
    Abar=Abar+Acells{k};
end
Abar=Abar/nC;

if max(real(eig(Abar)))<0
    P=rg_solve_lyap(Abar,eye(n));
    P=0.5*(P+P.');
    if min(eig(P))>0
        Pcand{end+1}=P; %#ok<AGROW>
    end
end

for k=1:nC
    if alphaSpec(k)<0
        P=rg_solve_lyap(Acells{k},eye(n));
        P=0.5*(P+P.');
        if min(eig(P))>0
            Pcand{end+1}=P; %#ok<AGROW>
        end
    end
end

% Always include a scaled identity as a fallback candidate.
Pcand{end+1}=eye(n);

if isempty(Pcand)
    error('viii_common_group_storage:NoCandidates','No SPD storage candidates.');
end

% Normalize candidate traces.
for k=1:numel(Pcand)
    Pcand{k}=Pcand{k}/(trace(Pcand{k})/n);
end

bestP=Pcand{1};
bestMu=local_worst_rate(bestP,Acells);

for k=2:numel(Pcand)
    mu=local_worst_rate(Pcand{k},Acells);
    if mu<bestMu
        bestMu=mu;
        bestP=Pcand{k};
    end
end

% Arithmetic mean candidate.
Pmean=zeros(n);
for k=1:numel(Pcand)
    Pmean=Pmean+Pcand{k};
end
Pmean=Pmean/numel(Pcand);
mu=local_worst_rate(Pmean,Acells);
if mu<bestMu
    bestMu=mu;
    bestP=Pmean;
end

% Deterministic pairwise convex-blending refinement.
for pass=1:P8.macro.maxBlendPasses
    improved=false;
    P0=bestP;
    mu0=bestMu;

    for k=1:numel(Pcand)
        for ia=1:numel(P8.macro.blendGrid)
            a=P8.macro.blendGrid(ia);
            Ptry=(1-a)*P0+a*Pcand{k};
            Ptry=0.5*(Ptry+Ptry.');
            Ptry=Ptry/(trace(Ptry)/n);

            mutry=local_worst_rate(Ptry,Acells);
            if mutry<bestMu-1e-10
                bestMu=mutry;
                bestP=Ptry;
                improved=true;
            end
        end
    end

    if ~improved || bestMu>=mu0-1e-10
        break;
    end
end

O.P=bestP;
O.mu=bestMu;
O.decayRaw=-bestMu;
O.maxSelfSpectralAbscissa=max(alphaSpec);
O.minEigP=min(eig(bestP));
O.condP=cond(bestP);
O.nCandidates=numel(Pcand);
end

function mu=local_worst_rate(P,Acells)
mu=-inf;
for kk=1:numel(Acells)
    mu=max(mu,viii_gen_rate(Acells{kk},P));
end
end
