function R = rotvec_to_dcm(rv)
%ROTVEC_TO_DCM Rodrigues map exp(hat(rv)); base-MATLAB implementation.
rv=rv(:); th=norm(rv);
if th < 1e-12
    K=hat3(rv);
    R=eye(3)+K+0.5*K*K;
else
    a=rv/th; K=hat3(a);
    R=eye(3)+sin(th)*K+(1-cos(th))*(K*K);
end
end

function H=hat3(a)
H=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
