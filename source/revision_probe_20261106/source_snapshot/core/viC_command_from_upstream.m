function O=viC_command_from_upstream(z,R,A)
%VIC_COMMAND_FROM_UPSTREAM Command attitude from Section V/VI-A upstream state.
% z=[xf(6); y(3)].

xf=z(1:6);
y=z(7:9);
C=viA_feedback_command(y,xf,R,A);

O.C=C;
O.Cc=viB_command_attitude(y(2),y(3),C.mu,C.alpha);
O.gamma=y(2);
O.chi=y(3);
O.mu=C.mu;
O.alpha=C.alpha;
end
