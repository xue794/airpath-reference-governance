function eul = dcm_to_euler321(C)
%DCM_TO_EULER321 Local 3-2-1 Euler angles [phi theta psi]' from C_IB.
s=-C(3,1); s=max(-1,min(1,s));
theta=asin(s);
if abs(cos(theta)) < 1e-10
    error('dcm_to_euler321:Singularity','Euler extraction near pitch singularity.');
end
phi=atan2(C(3,2),C(3,3));
psi=atan2(C(2,1),C(1,1));
eul=[phi;theta;psi];
end
