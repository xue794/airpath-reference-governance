function z=cl6_wind_profile(t,windMode,C)
%CL6_WIND_PROFILE Inertial NED wind for the integrated 6DOF cases.
% Output z=[w;wdot].

switch lower(windMode)
    case 'still'
        w=zeros(3,1);
        wd=zeros(3,1);

    case 'constant'
        w=C.E.w0;
        wd=zeros(3,1);

    case 'varying'
        [w,wd]=wind_estimator_test_wind(t,'varying');

    case 'certified_small'
        w=[4;-2;0];
        wd=zeros(3,1);

    case 'monte_carlo'
        if ~isfield(C,'MC') || ~isfield(C.MC,'wind')
            error('cl6_wind_profile:MC','Monte-Carlo wind parameters are missing.');
        end
        W=C.MC.wind;
        arg=W.omega(:)*t+W.phase(:);
        w=W.mean(:)+W.amp(:).*sin(arg);
        wd=W.amp(:).*W.omega(:).*cos(arg);

    otherwise
        error('cl6_wind_profile:Mode','Unknown wind mode: %s',windMode);
end

z=[w;wd];
end
