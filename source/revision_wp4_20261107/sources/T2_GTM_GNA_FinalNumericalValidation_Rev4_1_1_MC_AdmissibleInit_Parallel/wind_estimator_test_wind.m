function [w,wd] = wind_estimator_test_wind(t,mode)
%WIND_ESTIMATOR_TEST_WIND Deterministic wind profiles used in WE-P tests.
if nargin<2, mode='varying'; end
switch lower(mode)
    case 'constant'
        w=[8;-4;1]; wd=zeros(3,1);
    case 'varying'
        w=[6+2*sin(0.25*t); ...
           -3+1.5*sin(0.40*t+0.5); ...
           0.8*sin(0.20*t)];
        wd=[0.5*cos(0.25*t); ...
            0.6*cos(0.40*t+0.5); ...
            0.16*cos(0.20*t)];
    otherwise
        error('wind_estimator_test_wind:Mode','Unknown mode.');
end
end
