function mu=viii_gen_rate(A,P)
%VIII_GEN_RATE Largest generalized symmetric derivative rate.
%
% H=0.5(A'P+PA).  If mu<0 then
%   H <= mu P
% and V=0.5*x'Px has root-energy decay at least -mu.

H=0.5*(A.'*P+P*A);
H=0.5*(H+H.');
P=0.5*(P+P.');
ev=eig(H,P);
mu=max(real(ev));
end
