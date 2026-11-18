function G=viB_geometry(gamma,chi,mu,alphaC,Qf,QR)
%VIB_GEOMETRY Complete Section VI-B relative-attitude geometry.
%
% Qf = Cc'*Cr
% QR = Cr'*C
% Qc = Qf*QR = Cc'*C

B=viA_airpath_basis(gamma,chi);
Cc=viB_command_attitude(gamma,chi,mu,alphaC);

Cr=Cc*Qf;
C=Cr*QR;
Qc=Qf*QR;

bc=Cc.'*B.t;
A=viB_airdata_from_geometry(C,B.t);

G.B=B;
G.Cc=Cc; G.Cr=Cr; G.C=C;
G.Qf=Qf; G.QR=QR; G.Qc=Qc;
G.bc=bc; G.b=A.b;
G.alpha=A.alpha; G.beta=A.beta; G.kappa=A.kappa;
G.thetaF=dcm_distance(Cc,Cr);
G.thetaR=dcm_distance(Cr,C);
G.thetaC=dcm_distance(Cc,C);
end
