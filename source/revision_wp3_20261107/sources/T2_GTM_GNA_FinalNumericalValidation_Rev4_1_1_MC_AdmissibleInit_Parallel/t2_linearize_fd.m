function A = t2_linearize_fd(fun,x0)
%T2_LINEARIZE_FD Central finite-difference state Jacobian.
n=numel(x0); A=zeros(n);
for i=1:n
    hi=1e-6*max(1,abs(x0(i)));
    xp=x0; xm=x0; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi;
    A(:,i)=(fun(xp)-fun(xm))/(2*hi);
end
end
