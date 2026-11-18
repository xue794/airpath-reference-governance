function [A,coef] = fit_sinusoid_amplitude(t,y,Omega)
%FIT_SINUSOID_AMPLITUDE Least-squares steady-state sinusoid amplitude.
% Fits y = c0 + cs*sin(Omega*t)+cc*cos(Omega*t).
t=t(:); y=y(:);
X=[ones(size(t)),sin(Omega*t),cos(Omega*t)];
coef=X\y;
A=hypot(coef(2),coef(3));
end
