function M=viD_rate_moment_decomp(V,alpha,R)
%VID_RATE_MOMENT_DECOMP T-2 GNA rate-moment split K=-D+S.
%
% The exact GNA rate-dependent aerodynamic moment is
%    tau_rate = Komega * omega
% and this module decomposes it as
%    Komega = -Domega + Somega,
% where Domega=Domega' and Somega=-Somega'.
%
% The symmetric dissipative part is the rotational PHS damping port.
% The skew part is NOT inserted into the nominal PHS interconnection; it is
% compensated explicitly before allocation in Section VII.

P=R.P;
th=P.theta;
qS=0.5*R.rho*P.S*V^2;

K=zeros(3,3);

% Roll / yaw rate derivatives.
K(1,1)=qS*P.b*th(25)*P.b/(2*V);
K(1,3)=qS*P.b*th(26)*P.b/(2*V);
K(3,1)=qS*P.b*th(40)*P.b/(2*V);
K(3,3)=qS*P.b*th(41)*P.b/(2*V);

% Pitch-rate derivative, including the alpha-dependent GNA terms.
Cmq=th(31)+th(33)*alpha+th(34)*alpha^2+th(36)*alpha^3;
K(2,2)=qS*P.cbar*Cmq*P.cbar/(2*V);

Domega=-0.5*(K+K.');
Somega= 0.5*(K-K.');

M.K=K;
M.D=Domega;
M.S=Somega;
M.minEigD=min(eig(Domega));
M.maxEigD=max(eig(Domega));
end
