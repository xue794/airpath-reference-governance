function B = viA_airpath_basis(gamma,chi)
%VIA_AIRPATH_BASIS Air-relative tangent/climb/course basis in inertial NED.

B.t = [cos(gamma)*cos(chi);
       cos(gamma)*sin(chi);
      -sin(gamma)];

B.nGamma = [-sin(gamma)*cos(chi);
            -sin(gamma)*sin(chi);
            -cos(gamma)];

B.nChi = [-sin(chi);
           cos(chi);
           0];
end
