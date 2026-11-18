function C = quat_to_dcm(q)
%QUAT_TO_DCM Body-to-inertial DCM from scalar-first quaternion.
% q represents the active rotation B -> I and is normalized internally.
q=quat_normalize(q);
w=q(1); x=q(2); y=q(3); z=q(4);
C=[1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w); ...
   2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w); ...
   2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)];
end
