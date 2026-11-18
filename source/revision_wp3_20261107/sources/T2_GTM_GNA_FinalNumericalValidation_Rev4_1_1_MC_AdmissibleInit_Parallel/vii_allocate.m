function O=vii_allocate(tauC,V,alpha,beta,omega,R,P7)
%VII_ALLOCATE Convert virtual PHS torque to physical surface command.
%
% Physical T-2 GNA moment:
%   tau_a = tau0 - D omega + S omega + B delta.
%
% VI-D nominal target:
%   -D omega + tau_c.
%
% Therefore
%   B delta = tau_req
%   tau_req = tau_c - tau0 - S omega - tau_prop_nom.
%
% The present study plant uses point thrust through the CG, so tau_prop_nom=0.

M=vii_moment_components(V,alpha,beta,R);

tauPropNom=zeros(3,1);
tauReq=tauC-M.tau0-M.S*omega-tauPropNom;

A=vii_box_allocator(M.B,tauReq,P7.deltaMin,P7.deltaMax,P7.Wtau,P7.allocTol);

O.tauReq=tauReq;
O.deltaCmd=A.delta;
O.dAlloc=A.residual;
O.cost=A.cost;
O.active=A.active;
O.status=A.status;
O.M=M;
O.tauPropNom=tauPropNom;
end
