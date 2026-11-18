function C = euler321_dcm(phi,theta,psi)
%EULER321_DCM Body-to-NED-inertial DCM C_IB = Rz(psi) Ry(theta) Rx(phi).
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta); cy=cos(psi); sy=sin(psi);
C = [cy*ct, cy*st*sp-sy*cp, cy*st*cp+sy*sp; ...
     sy*ct, sy*st*sp+cy*cp, sy*st*cp-cy*sp; ...
     -st,   ct*sp,          ct*cp];
end
