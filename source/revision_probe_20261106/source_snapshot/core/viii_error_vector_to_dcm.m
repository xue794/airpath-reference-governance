function Q=viii_error_vector_to_dcm(e)
%VIII_ERROR_VECTOR_TO_DCM Local inverse of e=0.5(Q-Q')^vee.
%
% On the certified attitude tubes theta<pi/2,
%   e = sin(theta)*axis
% has the unique principal inverse used here.

n=norm(e);
if n<1e-14
    Q=eye(3);
    return;
end
if n>=1
    error('viii_error_vector_to_dcm:Domain', ...
        'Local SO(3) error vector norm must be < 1.');
end
theta=asin(n);
Q=rotvec_to_dcm((theta/n)*e);
end
