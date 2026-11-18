function R=viB_Rz(a)
%VIB_RZ Rotation consistent with euler321_dcm.
c=cos(a); s=sin(a);
R=[c -s 0; s c 0; 0 0 1];
end
