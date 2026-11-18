function Q = rg_matched_qp(x,u0,R,O,C)
%RG_MATCHED_QP Same 18 reference rows; weighted least-modification QP.
% Q=rg_matched_qp(x,u0,R) builds C via the frozen rg_constraint_affine.
% Fifth input C is for same-state/testing reuse. Its correspondence to x,R
% is the caller's responsibility. NEVER replace rg_scalar_governor.m in place.
% No Q.lambda or Q.lambdaJ is fabricated. W2-C must update diagnostic consumers.
if nargin<4 || isempty(O), O=w2b_options(); end
assert(abs(R.rhoB-O.RhoB)<=1e-12,'W2B:Rho','Expected frozen candidate rhoB=0.10.');
if nargin<5 || isempty(C), C=rg_constraint_affine(x,R); end
expected=w2b_row_names();
assert(isequal(C.names(:),expected(:)) && isequal(size(C.G),[18 3]), ...
    'W2B:Rows','Constraint row identity/order differs from the locked protocol.');
uB=rg_backup_control(x,R);
OO=O; OO.FeasibleStart=uB;
[u,diagQP]=w2b_qp3(u0,-C.G,C.c,O.Scale,OO);
Q=struct('Method','Matched-QP','Success',diagQP.Success,'u',u,'u0',u0(:), ...
    'uB',uB,'C',C,'phi',nan(18,1),'phi0',C.c+C.G*u0(:), ...
    'sigma',C.c+C.G*uB,'Objective',diagQP.Objective, ...
    'CorrectionNorm',NaN,'Intervened',false,'Solver',diagQP);
if ~diagQP.Success
    error('W2B:QPSolver','%s; no backup substitution. min backup row=%.17g', ...
        diagQP.Status,min(Q.sigma));
end
Q.phi=C.c+C.G*u;
Q.CorrectionNorm=norm((u-u0(:))./O.Scale);
Q.Intervened=Q.CorrectionNorm>1e-10;
end
