function [yw,nw] = wind_pseudomeasurement(vgI_m,CIB_m,z_m,w_true)
%WIND_PSEUDOMEASUREMENT Reconstruct inertial wind from measured ground velocity,
% attitude, and airdata.
%
%   y_w = v_g,m^I - C_IB,m * T(V_m,alpha_m,beta_m)
%
% If w_true is supplied, nw = yw-w_true is also returned.
if numel(z_m)~=3
    error('wind_pseudomeasurement:Size','z_m must be [V;alpha;beta].');
end
vaB_m = airdata_to_vaB(z_m(1),z_m(2),z_m(3));
yw = vgI_m(:) - CIB_m*vaB_m;
if nargin >= 4 && ~isempty(w_true)
    nw = yw - w_true(:);
else
    nw = [];
end
end
