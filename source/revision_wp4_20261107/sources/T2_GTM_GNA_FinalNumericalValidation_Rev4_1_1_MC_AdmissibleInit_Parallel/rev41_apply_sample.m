function C=rev41_apply_sample(C0,S,R)
%REV41_APPLY_SAMPLE Apply one randomized environment/plant realization.
%
% The controller and its nominal internal model C.R are deliberately NOT
% changed.  Only C.P (true plant), C.P7 actuator time constants, C.E.MC
% measurement waveform, C.MC wind, and initial-condition perturbations vary.

C=C0;

% True rigid-body plant.
C.P.m=C0.P.m*S.massScale;
C.P.W_N=C.P.m*C.P.g;

C.P.Ixx=C0.P.Ixx*S.inertiaScale(1);
C.P.Iyy=C0.P.Iyy*S.inertiaScale(2);
C.P.Izz=C0.P.Izz*S.inertiaScale(3);
C.P.Ixz=C0.P.Ixz*S.inertiaScale(4);
C.P.I=[C.P.Ixx 0 -C.P.Ixz; 0 C.P.Iyy 0; -C.P.Ixz 0 C.P.Izz];

if min(eig(0.5*(C.P.I+C.P.I.')))<=0
    error('Rev41:Inertia','Randomized inertia lost positive definiteness.');
end

C.P.theta=C0.P.theta.*S.aeroScale;

% True actuator time constants. Hard boxes/rate limits remain frozen.
C.P7.tauDelta=C0.P7.tauDelta.*S.servoTauScale;
C.P7.tauT=C0.P7.tauT*S.thrustTauScale;

% Smooth bounded measurement-error waveform.
C.E.MC.noiseScale=S.noiseScale;
C.E.MC.freqScale=S.noiseFreqScale;
C.E.MC.phase=S.noisePhase;

% Monte-Carlo environment and initial conditions.
C.MC.T=R.T;
C.MC.tSteady=R.tSteady;
C.MC.sampleIndex=S.index;
C.MC.initStrategy=R.initStrategy;
C.MC.wind.mean=S.windMean;
C.MC.wind.amp=S.windAmp;
C.MC.wind.omega=S.windOmega;
C.MC.wind.phase=S.windPhase;

C.MC.init.position=S.initPosition;
C.MC.init.vgB=S.initVgB;
C.MC.init.attitudeRotvec=S.initAttitudeRotvec;
C.MC.init.omega=S.initOmega;
C.MC.init.windEstimateError=S.initWhatError;
end
