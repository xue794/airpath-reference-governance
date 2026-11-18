function v = rg_airpath_to_velocity(c)
%RG_AIRPATH_TO_VELOCITY Convert [V;gamma;chi] to inertial NED air velocity.
V=c(1); gamma=c(2); chi=c(3);
v = V*[cos(gamma)*cos(chi);
       cos(gamma)*sin(chi);
      -sin(gamma)];
end
