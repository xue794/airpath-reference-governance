function tdot=viB_tangent_rate(gamma,chi,gammaDot,chiDot)
%VIB_TANGENT_RATE Time derivative of air-path tangent t(gamma,chi).
B=viA_airpath_basis(gamma,chi);
tdot=B.nGamma*gammaDot + cos(gamma)*B.nChi*chiDot;
end
