function W=viii_macro_small_gain_witness(Data,C,P8)
%VIII_MACRO_SMALL_GAIN_WITNESS Robust two-macro-subsystem diagnostic/witness.
%
% Every candidate partition is evaluated and retained in W.candidates.
% If no partition has positive sampled self decay, this function returns
% W.valid=false with diagnostics; it never throws "No macro partition evaluated."

J=viii_collect_macro_jacobians(Data,C,P8);
Parts=viii_macro_partitions();

Cand=struct('partition',{},'idx1',{},'idx2',{},'G1',{},'G2',{}, ...
            'a',{},'c',{},'rho',{},'selfValid',{},'smallGainValid',{});

bestIdx=1;
bestScore=inf;

for ip=1:numel(Parts)
    idx1=viii_blocks_to_indices(Parts(ip).g1Blocks);
    idx2=viii_blocks_to_indices(Parts(ip).g2Blocks);

    A11=cell(numel(J.A),1);
    A22=cell(numel(J.A),1);
    for k=1:numel(J.A)
        A11{k}=J.A{k}(idx1,idx1);
        A22{k}=J.A{k}(idx2,idx2);
    end

    G1=viii_common_group_storage(A11,P8);
    G2=viii_common_group_storage(A22,P8);

    selfValid=(G1.decayRaw>P8.macro.minDecay) && ...
              (G2.decayRaw>P8.macro.minDecay);

    if selfValid
        a1=P8.macro.selfDeflation*G1.decayRaw;
        a2=P8.macro.selfDeflation*G2.decayRaw;

        R1=chol(G1.P);
        R2=chol(G2.P);

        c12=0;
        c21=0;
        for k=1:numel(J.A)
            A12=J.A{k}(idx1,idx2);
            A21=J.A{k}(idx2,idx1);
            c12=max(c12,norm(R1*A12/R2,2));
            c21=max(c21,norm(R2*A21/R1,2));
        end

        c12=P8.macro.crossInflation*c12;
        c21=P8.macro.crossInflation*c21;

        D=diag([a1 a2]);
        G=[0 c12;c21 0];
        rho=max(abs(eig(D\G)));
    else
        a1=P8.macro.selfDeflation*G1.decayRaw;
        a2=P8.macro.selfDeflation*G2.decayRaw;
        c12=nan;
        c21=nan;
        rho=inf;
    end

    smallValid=selfValid && isfinite(rho) && ...
               rho<1-P8.macro.rhoMargin;

    Cand(ip).partition=Parts(ip);
    Cand(ip).idx1=idx1;
    Cand(ip).idx2=idx2;
    Cand(ip).G1=G1;
    Cand(ip).G2=G2;
    Cand(ip).a=[a1;a2];
    Cand(ip).c=[0 c12;c21 0];
    Cand(ip).rho=rho;
    Cand(ip).selfValid=selfValid;
    Cand(ip).smallGainValid=smallValid;

    % Ranking: valid small gain first, then finite rho, then best self-decay deficit.
    if smallValid
        score=rho;
    elseif isfinite(rho)
        score=10+rho;
    else
        deficit=max(0,P8.macro.minDecay-G1.decayRaw)+ ...
                max(0,P8.macro.minDecay-G2.decayRaw);
        score=1e6+deficit;
    end

    if ip==1 || score<bestScore
        bestScore=score;
        bestIdx=ip;
    end
end

best=Cand(bestIdx);

D=diag(best.a);
G=best.c;

if best.smallGainValid
    M=D-G;
    v=M\ones(2,1);
    w=M.'\ones(2,1);
    piw=w./v;
    Pi=diag(piw);
    Q=M.'*Pi+Pi*M;
    minEigQ=min(eig(0.5*(Q+Q.')));
else
    M=nan(2);
    v=nan(2,1);
    w=nan(2,1);
    piw=nan(2,1);
    Pi=nan(2);
    Q=nan(2);
    minEigQ=nan;
end

B=viii_blocks();
Pfull=zeros(B.n);
Pfull(best.idx1,best.idx1)=best.G1.P;
Pfull(best.idx2,best.idx2)=best.G2.P;

W=best;
W.valid=best.smallGainValid;
W.D=D;
W.G=G;
W.M=M;
W.v=v;
W.w=w;
W.pi=piw;
W.Pi=Pi;
W.Q=Q;
W.minEigQ=minEigQ;
W.Pfull=Pfull;
W.J=J;
W.nWitnessPoints=numel(J.A);
W.candidates=Cand;
W.bestIndex=bestIdx;
end
