function R=viB_Rx(a)
%VIB_RX Active/body-frame x-axis rotation used in C_IB products.
c=cos(a); s=sin(a);
R=[1 0 0; 0 c -s; 0 s c];
end
