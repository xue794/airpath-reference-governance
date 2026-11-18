function z0 = viA_initial_state(mode,R,A)
%VIA_INITIAL_STATE Combined governor + ideal air-path initial condition.
if strcmpi(mode,'wrap')
    [~,c0]=rg_raw_test_command(0,'wrap',R);
    xf0=[R.Vb;R.gammaB;c0(3);0;0;0];
    e0=A.e0Wrap;
else
    xf0=[R.Vb;R.gammaB;0;0;0;0];
    e0=A.e0;
end

y0=[xf0(1)+e0(1);
    xf0(2)+e0(2);
    xf0(3)+e0(3)];

z0=[xf0;y0];
end
