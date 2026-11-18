function O=cl6_allocate_fast(tauC,V,alpha,beta,omega,C)
%CL6_ALLOCATE_FAST Same Section-VII box-LS allocator with a feasible direct path.
%
% Most complete-loop commands are far inside the box.  A direct solve is
% therefore used first; the exact 27-active-set allocator is called only
% when the unconstrained solution leaves the box.

P7=C.P7;
M=vii_moment_components(V,alpha,beta,C.R);

tauPropNom=zeros(3,1);
tauReq=tauC-M.tau0-M.S*omega-tauPropNom;

d0=M.B\tauReq;

if all(d0>=P7.deltaMin-P7.allocTol) && all(d0<=P7.deltaMax+P7.allocTol)
    deltaCmd=min(max(d0,P7.deltaMin),P7.deltaMax);
    dAlloc=M.B*deltaCmd-tauReq;

    O.deltaCmd=deltaCmd;
    O.dAlloc=dAlloc;
    O.cost=norm(P7.Wtau*dAlloc)^2;
    O.active=false(3,1);
    O.status=zeros(3,1);
else
    A=vii_box_allocator(M.B,tauReq,P7.deltaMin,P7.deltaMax,P7.Wtau,P7.allocTol);
    O.deltaCmd=A.delta;
    O.dAlloc=A.residual;
    O.cost=A.cost;
    O.active=A.active;
    O.status=A.status;
end

O.tauReq=tauReq;
O.M=M;
O.tauPropNom=tauPropNom;
end
