function [vc,cTrue] = rg_raw_test_command(t,mode,R)
%RG_RAW_TEST_COMMAND Smooth raw air-relative velocity commands for Section V tests.
switch lower(mode)
    case 'small'
        V = R.Vb + 0.30*sin(0.05*t);
        gamma = deg2rad(0.50)*sin(0.04*t);
        chi = deg2rad(3.0)*sin(0.03*t);
    case 'aggressive'
        V = R.Vb + 3.0*sin(0.10*t);
        gamma = deg2rad(6.0)*sin(0.08*t);
        chi = deg2rad(40.0)*sin(0.07*t);
    case 'wrap'
        V = R.Vb;
        gamma = 0;
        chi = deg2rad(170.0) + deg2rad(0.50)*t;
    otherwise
        error('Unknown reference-governor command mode: %s',mode);
end
cTrue=[V;gamma;chi];
vc=rg_airpath_to_velocity(cTrue);
end
