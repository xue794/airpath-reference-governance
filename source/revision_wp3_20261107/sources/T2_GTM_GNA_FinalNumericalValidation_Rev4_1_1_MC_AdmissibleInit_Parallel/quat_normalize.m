function qn = quat_normalize(q)
%QUAT_NORMALIZE Normalize a scalar-first quaternion q=[q0 q1 q2 q3]'.
q=q(:); n=norm(q);
if n < 1e-12
    error('quat_normalize:ZeroNorm','Quaternion norm is too small.');
end
qn=q/n;
end
