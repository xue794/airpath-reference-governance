function y=w2e_quantiles(x,p)
x=x(isfinite(x));x=sort(x(:));n=numel(x);y=nan(size(p));if n==0,return;end
for k=1:numel(p),h=1+(n-1)*p(k);lo=floor(h);hi=ceil(h);if lo==hi,y(k)=x(lo);else,y(k)=x(lo)+(h-lo)*(x(hi)-x(lo));end;end
end
