function G=w2c_governor(x,u0,R,method,O,Ccon)
%W2C_GOVERNOR Select the frozen Proposed scalar method or Matched-QP.
% No QP lambda is fabricated. Both methods use the same 18-row constraint map.
if nargin<5 || isempty(O),O=w2c_options();end
if nargin<6,Ccon=[];end
switch lower(method)
    case 'proposed'
        G=rg_scalar_governor(x,u0,R);
        if isempty(Ccon),Ccon=G.C;end
        G.Method='Proposed';
        G.Success=true;
        G.phi0=Ccon.c+Ccon.G*u0(:);
        G.Objective=0.5*sum(((G.u-u0(:))./O.QP.Scale).^2);
        G.CorrectionNorm=norm((G.u-u0(:))./O.QP.Scale);
        G.Intervened=G.CorrectionNorm>O.InterventionTol;
        G.Solver=[];
    case 'matchedqp'
        G=rg_matched_qp(x,u0,R,O.QP,Ccon);
        G.Method='MatchedQP';
    otherwise
        error('W2C:Method','Unknown method: %s',method);
end
end
