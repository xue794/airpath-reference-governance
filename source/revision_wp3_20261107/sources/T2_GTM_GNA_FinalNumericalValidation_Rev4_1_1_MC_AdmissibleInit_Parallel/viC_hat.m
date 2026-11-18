function H=viC_hat(v)
%VIC_HAT 3-vector to skew-symmetric matrix.
v=v(:);
H=[0 -v(3) v(2);
   v(3) 0 -v(1);
  -v(2) v(1) 0];
end
