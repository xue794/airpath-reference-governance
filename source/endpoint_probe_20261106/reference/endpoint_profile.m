function cc=endpoint_profile(t,S,R)
% Common exploratory C2 two-stage request, from the legal trim.
a=blend((t-3)/8);b=blend((t-16)/3);
cc=[R.Vb+(S.speed-R.Vb)*a;S.gamma*pi/180*b;S.course*pi/180*b];
end
function a=blend(q)
q=min(1,max(0,q));a=q^3*(10+q*(-15+6*q));
end
