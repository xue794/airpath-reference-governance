function [s,I]=w2d_smooth_pulse(t,t0,Tr,Th)
%W2D_SMOOTH_PULSE C2 rise-hold-fall pulse and its time integral.
% Rise/fall use f(q)=10q^3-15q^4+6q^5. The final blend returns to zero.
if Tr<=0 || Th<0
    error('W2D:Pulse','Ramp duration must be positive and hold duration nonnegative.');
end
t1=t0+Tr;
t2=t1+Th;
t3=t2+Tr;
if t<=t0
    s=0;I=0;
elseif t<t1
    q=(t-t0)/Tr;
    f=10*q^3-15*q^4+6*q^5;
    F=2.5*q^4-3*q^5+q^6;
    s=f;I=Tr*F;
elseif t<=t2
    s=1;I=0.5*Tr+(t-t1);
elseif t<t3
    q=(t-t2)/Tr;
    f=10*q^3-15*q^4+6*q^5;
    F=2.5*q^4-3*q^5+q^6;
    s=1-f;I=0.5*Tr+Th+Tr*(q-F);
else
    s=0;I=Tr+Th;
end
end
