function q = w2f_percentile(x,p)
%W2F_PERCENTILE Linear interpolation at h=1+(n-1)p; base MATLAB only.
x=x(:);x=x(isfinite(x));
assert(~isempty(x),'W2F:Stats','No finite observations.');
assert(p>=0 && p<=1,'W2F:Stats','p must be in [0,1].');
x=sort(x);n=numel(x);h=1+(n-1)*p;k=floor(h);d=h-k;
if k>=n,q=x(n);else,q=(1-d)*x(k)+d*x(k+1);end
end
