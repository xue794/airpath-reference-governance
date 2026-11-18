function [rd,rdDot] = smc_guidance_reference(t)
%SMC_GUIDANCE_REFERENCE Smooth NED ground-trajectory test reference.
% NED coordinates: D is positive down.
%
% This trajectory is used only to exercise the guidance implementation.
% The error dynamics are reference-invariant after the exact cancellations.

rd = [20.0*t;
      30.0*sin(0.05*t);
     -365.76 - 8.0*sin(0.04*t)];

rdDot = [20.0;
          1.50*cos(0.05*t);
         -0.32*cos(0.04*t)];
end
