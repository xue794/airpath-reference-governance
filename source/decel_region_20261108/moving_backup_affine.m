function [C,uBackup,h,center] = moving_backup_affine(t,xf,R)
% Fixed 22 s/0.01 certified candidate. xf = [V,gamma,chi,nuV,nuGamma,nuChi].
% Pair the returned affine rows with the returned feasible backup input.
if t < 0 || t > 55
    error('MovingBackup:time','Certified implementation domain is 0..55 s.');
end
if t<=3
    v=R.Vb; n=0; a=0;
elseif t<=25
    [v,n,a]=blend(R.Vb,31.25,(t-3)/22,22);
elseif t<=35
    v=31.25; n=0; a=0;
else
    [v,n,a]=blend(31.25,R.Vb,(t-35)/20,20);
end
center=[v;0;n;0;0]; z=xf([1 2 4 5 6]);e=z-center;
y=R.PB*e; h=.01-e'*y;
uBackup=[a-R.kVp*e(1)-R.kVd*e(3);
        -R.kGammaP*e(2)-R.kGammaD*e(4);
        -R.kChiD*e(5)];
C=rg_constraint_affine(xf,R);
C.c(18)=R.kB*h-2*y(1)*e(3)-2*y(2)*e(4)+2*y(3)*a;
C.G(18,:)=-2*y(3:5)';
C.hB=h;
end
function [v,n,a]=blend(v0,v1,q,d)
v=v0+(v1-v0)*q^3*(10-15*q+6*q^2);
n=(v1-v0)*30*q^2*(1-q)^2/d;
a=(v1-v0)*60*q*(1-q)*(1-2*q)/d^2;
end
