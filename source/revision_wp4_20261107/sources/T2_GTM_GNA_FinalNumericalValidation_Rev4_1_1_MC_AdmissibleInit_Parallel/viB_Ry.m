function R=viB_Ry(a)
%VIB_RY Rotation consistent with euler321_dcm: Rz(psi)Ry(theta)Rx(phi).
c=cos(a); s=sin(a);
R=[c 0 s; 0 1 0; -s 0 c];
end
