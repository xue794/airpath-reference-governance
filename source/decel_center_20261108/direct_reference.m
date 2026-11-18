function [xf,accel]=direct_reference(t,R)
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
xf=[v;0;0;n;0;0];accel=[a;0;0];
end
function [v,n,a]=blend(v0,v1,q,d)
v=v0+(v1-v0)*q^3*(10-15*q+6*q^2);
n=(v1-v0)*30*q^2*(1-q)^2/d;
a=(v1-v0)*60*q*(1-q)*(1-2*q)/d^2;
end
