function S = viA_case_signal(t,mode,E)
%VIA_CASE_SIGNAL Reference-governor mode and wind signal for VI-A tests.
switch lower(mode)
    case 'still'
        S.refMode='small';
        S.w=zeros(3,1); S.wdot=zeros(3,1);

    case 'constant'
        S.refMode='small';
        S.w=E.w0; S.wdot=zeros(3,1);

    case 'varying'
        S.refMode='small';
        [S.w,S.wdot]=wind_estimator_test_wind(t,'varying');

    case 'aggressive'
        % Feasibility-transfer test isolates reference motion from wind acceleration.
        S.refMode='aggressive';
        S.w=zeros(3,1); S.wdot=zeros(3,1);

    case 'wrap'
        S.refMode='wrap';
        S.w=zeros(3,1); S.wdot=zeros(3,1);

    otherwise
        error('Unknown VI-A mode: %s',mode);
end
end
