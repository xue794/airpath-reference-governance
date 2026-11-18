function A=d0_air(t,xp,windFun)
z=windFun(t);
wI=z(1:3);
wDotI=z(4:6);

CIB=quat_to_dcm(xp(7:10));
vgB=xp(4:6);
vgI=CIB*vgB;

vaB=vgB-CIB.'*wI;
vaI=CIB*vaB;

V=norm(vaI);
if V<1e-6
    error('cl6_airpath_state:V','Airspeed too small.');
end

Vh=hypot(vaI(1),vaI(2));
gamma=atan2(-vaI(3),Vh);
chi=atan2(vaI(2),vaI(1));

alpha=atan2(vaB(3),vaB(1));
beta=asin(max(-1,min(1,vaB(2)/V)));

A.CIB=CIB;
A.vgB=vgB;
A.vgI=vgI;
A.vaB=vaB;
A.vaI=vaI;
A.V=V;
A.gamma=gamma;
A.chi=chi;
A.alpha=alpha;
A.beta=beta;
A.kappa=cos(alpha)*cos(beta);
A.wI=wI;
A.wDotI=wDotI;
A.omega=xp(11:13);
A.h=-xp(3);
end
