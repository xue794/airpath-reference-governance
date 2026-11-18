function q = euler321_to_quat(phi,theta,psi)
%EULER321_TO_QUAT Scalar-first quaternion for C_IB=Rz(psi)Ry(theta)Rx(phi).
hp=0.5*phi; ht=0.5*theta; hs=0.5*psi;
cp=cos(hp); sp=sin(hp); ct=cos(ht); st=sin(ht); cs=cos(hs); ss=sin(hs);
q=[cs*ct*cp + ss*st*sp; ...
   cs*ct*sp - ss*st*cp; ...
   cs*st*cp + ss*ct*sp; ...
   ss*ct*cp - cs*st*sp];
q=quat_normalize(q);
end
