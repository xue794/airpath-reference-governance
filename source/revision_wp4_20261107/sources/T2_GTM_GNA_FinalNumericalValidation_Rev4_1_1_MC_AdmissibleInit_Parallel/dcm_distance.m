function a = dcm_distance(C1,C2)
%DCM_DISTANCE Principal rotation angle between two DCMs [rad].
Q=C1.'*C2;
c=(trace(Q)-1)/2; c=max(-1,min(1,c));
a=acos(c);
end
