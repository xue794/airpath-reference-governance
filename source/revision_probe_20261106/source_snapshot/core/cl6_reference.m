function [rd,rdDot]=cl6_reference(t,refMode,C)
%CL6_REFERENCE Smooth NED ground-position references for full 6DOF tests.
% NED: D positive down.

h0=C.R.h_ref;

switch lower(refMode)
    case 'hold'
        rd=[C.VgHold*t;
            0;
           -h0];
        rdDot=[C.VgHold;0;0];

    case 'constant'
        rd=[C.VgConstant*t;
            40*sin(0.01*t);
           -h0];
        rdDot=[C.VgConstant;
               0.40*cos(0.01*t);
               0];

    case 'varying'
        rd=[C.VgVarying*t;
            50*sin(0.01*t);
           -h0-6*sin(0.008*t)];
        rdDot=[C.VgVarying;
               0.50*cos(0.01*t);
              -0.048*cos(0.008*t)];

    case 'certified_small'
        % Strictly interior mild maneuver around the T-2 trim airspeed.
        % The associated constant wind is [4;-2;0] m/s, so the mean
        % air-relative velocity is [V_b;0;0].  The perturbations are
        % intentionally small enough to test the LOCAL theorem, not the
        % aggressive robustness envelope.
        w0=[4;-2;0];

        % Numerical-validation-only trajectory amplitude scale.
        % This field never changes controller/plant/actuator parameters.
        if isfield(C,'CertScale')
            sc=C.CertScale;
        else
            sc=1.0;
        end

        aN=sc*0.20; omN=0.025;
        aE=sc*0.60; omE=0.020;
        aD=sc*0.08; omD=0.018;

        rd=[ (C.R.Vb+w0(1))*t + (aN/omN)*(1-cos(omN*t));
             w0(2)*t          + (aE/omE)*(1-cos(omE*t));
            -h0               + (aD/omD)*(1-cos(omD*t)) ];

        rdDot=[ C.R.Vb+w0(1) + aN*sin(omN*t);
                w0(2)       + aE*sin(omE*t);
                               aD*sin(omD*t) ];

    case 'wrap'
        sp=C.wrapSpeed;
        chi0=C.wrapChi0;
        om=C.wrapChiRate;
        chi=chi0+om*t;

        rd=[ sp/om*(sin(chi)-sin(chi0));
             sp/om*(-cos(chi)+cos(chi0));
            -h0 ];
        rdDot=[sp*cos(chi);
               sp*sin(chi);
               0];

    otherwise
        error('cl6_reference:Mode','Unknown reference mode: %s',refMode);
end
end
