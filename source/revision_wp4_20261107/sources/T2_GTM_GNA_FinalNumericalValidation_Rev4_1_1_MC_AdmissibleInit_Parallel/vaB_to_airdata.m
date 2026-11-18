function z = vaB_to_airdata(vaB)
%VAB_TO_AIRDATA Convert body-frame air-relative velocity to [V;alpha;beta].
V = norm(vaB);
if V <= 0
    error('vaB_to_airdata:ZeroSpeed','Airspeed must be positive.');
end
alpha = atan2(vaB(3),vaB(1));
sarg = max(-1,min(1,vaB(2)/V));
beta = asin(sarg);
z = [V;alpha;beta];
end
