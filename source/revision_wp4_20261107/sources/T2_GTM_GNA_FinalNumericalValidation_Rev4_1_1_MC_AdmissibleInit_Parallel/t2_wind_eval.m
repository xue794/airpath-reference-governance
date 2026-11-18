function [wI,wDotI] = t2_wind_eval(t,U)
%T2_WIND_EVAL Uniform inertial wind evaluator.
% Backward-compatible forms:
%   U.wI = constant 3x1 wind [m/s]
% Optional time-varying form:
%   U.wind_fun = @(t) [wN;wE;wD;wdotN;wdotE;wdotD]
% If wind_fun returns only 3 entries, wDotI is set to zero.
if isfield(U,'wind_fun') && ~isempty(U.wind_fun)
    z=U.wind_fun(t); z=z(:);
    if numel(z)==3
        wI=z; wDotI=zeros(3,1);
    elseif numel(z)==6
        wI=z(1:3); wDotI=z(4:6);
    else
        error('t2_wind_eval:Size','U.wind_fun(t) must return 3 or 6 entries.');
    end
elseif isfield(U,'wI')
    wI=U.wI(:);
    if numel(wI)~=3, error('t2_wind_eval:Size','U.wI must have 3 entries.'); end
    if isfield(U,'wDotI')
        wDotI=U.wDotI(:);
    else
        wDotI=zeros(3,1);
    end
else
    wI=zeros(3,1); wDotI=zeros(3,1);
end
end
