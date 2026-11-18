function O = w2f_options()
%W2F_OPTIONS Frozen publication timing protocol for MATLAB R2018a.
% Do not tune these settings after viewing timing results.
B=w2b_options();
O=B;
O.Version='W2F_R2018a_v1.0';
O.Protocol='P1_Aerospace_W2F_PublicationTiming_rho010_v1_0';
O.RhoB=0.10;
O.Scale=[2;8*pi/180;12*pi/180];
O.WarmupPasses=3;
O.NumBlocks=12;
O.PassesPerBlock=5;
O.PermutationStride=137;   % gcd(137,425)=1; deterministic, no RNG use.
O.BlockShift=37;
O.RequireNoParallelPool=true;
O.ForceSingleComputeThread=true;
O.CorrectnessTol=1e-10;
O.ReferenceMatchTol=1e-8;
O.PauseBeforeTiming_s=1.0;
O.PublicationBenchmark=true;
end
