function z=confirm_wind(t)
q=min(1,max(0,t/3));s=q^3*(10-15*q+6*q^2);sd=10*q^2*(1-q)^2;
[w,wd]=wind_estimator_test_wind(t,'varying');
z=.1*[s*w;sd*w+s*wd];
end
