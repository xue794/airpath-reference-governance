function C=viD_damping_certificate(R,Dp)
%VID_DAMPING_CERTIFICATE Deterministic grid bounds for D_omega,a.

minEigD=inf;
maxNormD=0;
maxDecomp=0;
maxSkew=0;
maxSym=0;

Vgrid=linspace(R.Vmin,R.Vmax,Dp.scanNV);
agrid=linspace(R.alphaPhysicalMin,R.alphaPhysicalMax,Dp.scanNA);

for i=1:numel(Vgrid)
    for j=1:numel(agrid)
        M=viD_rate_moment_decomp(Vgrid(i),agrid(j),R);
        minEigD=min(minEigD,min(eig(M.D)));
        maxNormD=max(maxNormD,norm(M.D,2));
        maxDecomp=max(maxDecomp,norm(M.K+M.D-M.S,'fro'));
        maxSkew=max(maxSkew,norm(M.S+M.S.','fro'));
        maxSym=max(maxSym,norm(M.D-M.D.','fro'));
    end
end

C.minEigD=minEigD;
C.maxNormD=maxNormD;
C.maxDecomp=maxDecomp;
C.maxSkew=maxSkew;
C.maxSym=maxSym;
end
