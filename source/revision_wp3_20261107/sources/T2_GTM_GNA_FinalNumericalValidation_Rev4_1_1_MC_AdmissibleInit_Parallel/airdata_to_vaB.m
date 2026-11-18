function vaB = airdata_to_vaB(V,alpha,beta)
%AIRDATA_TO_VAB Reconstruct body-frame air-relative velocity from airdata.
% Inputs: V [m/s], alpha [rad], beta [rad].
vaB = V*[cos(alpha)*cos(beta); sin(beta); sin(alpha)*cos(beta)];
end
