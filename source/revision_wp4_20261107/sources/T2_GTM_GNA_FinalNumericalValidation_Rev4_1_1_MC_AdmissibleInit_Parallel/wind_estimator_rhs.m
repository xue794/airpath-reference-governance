function whatDot = wind_estimator_rhs(~,what,yw,Lw)
%WIND_ESTIMATOR_RHS First-order inertial wind estimator.
%   dot(hat w) = Lw (y_w - hat w)
whatDot = Lw*(yw(:)-what(:));
end
