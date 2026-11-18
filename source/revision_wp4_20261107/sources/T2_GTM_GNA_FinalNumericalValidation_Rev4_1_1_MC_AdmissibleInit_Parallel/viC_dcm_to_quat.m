function q=viC_dcm_to_quat(C)
%VIC_DCM_TO_QUAT Scalar-first quaternion for a body-to-inertial DCM.
% Base-MATLAB implementation, MATLAB R2018a compatible.

tr=trace(C);
if tr>0
    s=sqrt(tr+1.0)*2;
    qw=0.25*s;
    qx=(C(3,2)-C(2,3))/s;
    qy=(C(1,3)-C(3,1))/s;
    qz=(C(2,1)-C(1,2))/s;
elseif (C(1,1)>C(2,2)) && (C(1,1)>C(3,3))
    s=sqrt(1.0+C(1,1)-C(2,2)-C(3,3))*2;
    qw=(C(3,2)-C(2,3))/s;
    qx=0.25*s;
    qy=(C(1,2)+C(2,1))/s;
    qz=(C(1,3)+C(3,1))/s;
elseif C(2,2)>C(3,3)
    s=sqrt(1.0+C(2,2)-C(1,1)-C(3,3))*2;
    qw=(C(1,3)-C(3,1))/s;
    qx=(C(1,2)+C(2,1))/s;
    qy=0.25*s;
    qz=(C(2,3)+C(3,2))/s;
else
    s=sqrt(1.0+C(3,3)-C(1,1)-C(2,2))*2;
    qw=(C(2,1)-C(1,2))/s;
    qx=(C(1,3)+C(3,1))/s;
    qy=(C(2,3)+C(3,2))/s;
    qz=0.25*s;
end
q=quat_normalize([qw;qx;qy;qz]);
if q(1)<0
    q=-q;
end
end
