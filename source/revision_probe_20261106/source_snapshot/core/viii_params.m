function P8=viii_params()
%VIII_PARAMS Final composite Lyapunov / small-gain certificate parameters.
% MATLAB R2018a compatible; base MATLAB only.
%
% Rev.1.1.1 packaging cleanup. No controller/plant parameter is changed.

C=cl6_params();

P8.revision='Section VIII Final Composite Lyapunov / Small-Gain Certificate Rev.1.1.2';
P8.labels={'z','s','w','p','G','F','R','A'};
P8.nSub=8;
P8.Vstar=C.R.Vb;

P8.scale.z=[30;30;15];
P8.scale.s=[2;2;1];
P8.scale.w=[1;1;0.5];
P8.scale.p=[1;1;1];
P8.scale.G=[10;deg2rad(10);deg2rad(15); ...
            2;deg2rad(8);deg2rad(12)];
P8.scale.F=[sin(C.B.thetaFMax)*ones(3,1); ...
            0.08*ones(3,1)];
P8.scale.R=[sin(C.B.thetaRMax)*ones(3,1); ...
            0.05*ones(3,1)];
P8.scale.A=[deg2rad(2)*ones(3,1);2.0];

P8.scaleVector=[P8.scale.z;P8.scale.s;P8.scale.w;P8.scale.p; ...
                P8.scale.G;P8.scale.F;P8.scale.R;P8.scale.A];

P8.flowDt=2e-6;
P8.jacStep=2e-4;

% Historical eight-block diagnostics only.
P8.crossInflation=1.05;
P8.selfDeflation=0.95;
P8.forcingInflation=1.20;

% Feasible sampling.
P8.sample.hold=[5 20 40 55];
P8.sample.varying=[5 20 40 60 80 100 115];
P8.sample.noisy=[10 40 70 100 115];
P8.sample.actualVarying=[5 15 25 35 45 55 65 75 85 95 105 115];
P8.sample.actualNoisy=[10 20 30 40 50 60 70 80 90 100 110 115];

P8.roundTripTol=5e-6;
P8.forbiddenBlockTol=5e-4;
P8.storageRatioTol=5e-3;
P8.residualIdentityTol=2e-10;
P8.smallGainMargin=1e-3;
P8.diagStabilityTol=1e-9;

P8.kappaDelta=min(1./C.P7.tauDelta, ...
    C.P7.deltaRateMax./(C.P7.deltaMax-C.P7.deltaMin));
P8.kappaT=min(1/C.P7.tauT, ...
    C.P7.TRateMax/(C.P7.Tmax-C.P7.Tmin));
P8.kappaA=min([P8.kappaDelta;P8.kappaT]);
P8.Nsector=401;

% Rev.1.1 macro-composite witness.
P8.macro.selfDeflation=0.90;
P8.macro.crossInflation=1.10;
P8.macro.forcingInflation=1.20;
P8.macro.minDecay=1e-5;
P8.macro.rhoMargin=1e-3;
P8.macro.blendGrid=[0.05 0.10 0.20 0.35 0.50 0.65 0.80 0.90 0.95];
P8.macro.maxBlendPasses=12;


% Rev.1.1.2 robust certificate-harness settings.
P8.cache.enable=true;
P8.cache.file='viii_preflight_cache_rev112.mat';
P8.direct.selfDeflation=0.90;
P8.direct.minDecay=1e-5;
P8.direct.forcingInflation=1.20;

end
