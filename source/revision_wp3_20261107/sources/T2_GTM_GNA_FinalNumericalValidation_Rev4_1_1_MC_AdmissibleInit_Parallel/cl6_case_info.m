function K=cl6_case_info(mode,C)
%CL6_CASE_INFO Scenario definition for complete 6DOF integration.

switch lower(mode)
    case 'hold'
        K.refMode='hold';
        K.windMode='still';
        K.addNoise=false;
        K.coldWindStart=false;
        K.T=C.Thold;
        K.tSteady=C.tSteadyHold;

    case 'constant_cold'
        K.refMode='constant';
        K.windMode='constant';
        K.addNoise=false;
        K.coldWindStart=true;
        K.T=C.Tconstant;
        K.tSteady=C.tSteadyConstant;

    case 'varying'
        K.refMode='varying';
        K.windMode='varying';
        K.addNoise=false;
        K.coldWindStart=false;
        K.T=C.Tvarying;
        K.tSteady=C.tSteadyVarying;

    case 'noisy'
        K.refMode='varying';
        K.windMode='varying';
        K.addNoise=true;
        K.coldWindStart=false;
        K.T=C.Tnoisy;
        K.tSteady=C.tSteadyNoisy;

    case 'certified_small'
        % Mild nontrivial trajectory used ONLY to instantiate the local
        % Section-VIII theorem conditions. Controller parameters are unchanged.
        K.refMode='certified_small';
        K.windMode='certified_small';
        K.addNoise=false;
        K.coldWindStart=false;
        K.T=120.0;
        K.tSteady=30.0;

    case 'wrap'
        K.refMode='wrap';
        K.windMode='still';
        K.addNoise=false;
        K.coldWindStart=false;
        K.T=C.Twrap;
        K.tSteady=30.0;

    case 'monte_carlo'
        % Rev.4.1 extended randomized robustness case.
        % The controller remains nominal/frozen.  Only the plant,
        % environment, measurement-error waveform, actuator time constants,
        % and initial conditions are randomized through C.MC.
        if ~isfield(C,'MC')
            error('cl6_case_info:MC','mode=monte_carlo requires C.MC.');
        end
        K.refMode='varying';
        K.windMode='monte_carlo';
        K.addNoise=true;
        K.coldWindStart=false;
        K.T=C.MC.T;
        K.tSteady=C.MC.tSteady;

    otherwise
        error('cl6_case_info:Mode','Unknown complete-6DOF mode: %s',mode);
end

K.mode=mode;
end
