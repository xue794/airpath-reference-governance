function G=w2e_governor(x,u0,R,method,O,Ccon)
%W2E_GOVERNOR Frozen Proposed selector or frozen W2-B Matched-QP selector.
if nargin<5||isempty(O),O=w2e_options();end
if nargin<6,Ccon=[];end
switch lower(method)
 case 'proposed'
  G=rg_scalar_governor(x,u0,R);if isempty(Ccon),Ccon=G.C;end
  G.Method='Proposed';G.Success=true;G.phi0=Ccon.c+Ccon.G*u0(:);
  G.Objective=0.5*sum(((G.u-u0(:))./O.QP.Scale).^2);
  G.CorrectionNorm=norm((G.u-u0(:))./O.QP.Scale);
  G.Intervened=G.CorrectionNorm>O.InterventionTol;G.Solver=[];
 case 'matchedqp'
  G=rg_matched_qp(x,u0,R,O.QP,Ccon);G.Method='MatchedQP';
 otherwise,error('W2E:Method','Unknown method: %s',method);
end
end
