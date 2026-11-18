function S=simulate_complete_6dof_case(mode,C)
%SIMULATE_COMPLETE_6DOF_CASE Integrate one complete closed-loop scenario.

if nargin<2 || isempty(C)
    C=cl6_params();
end

K=cl6_case_info(mode,C);
[x0,I0]=cl6_initial_state(mode,C);

tgrid=(0:C.sample_dt:K.T).';
if tgrid(end)<K.T
    tgrid=[tgrid;K.T];
end

opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep);

[t,x]=ode45(@(tt,xx)cl6_rhs(tt,xx,mode,C),tgrid,x0,opts);

Y=cl6_outputs(t,x,mode,C);

S=struct('t',t,'x',x,'Y',Y,'mode',mode,'Case',K,'I0',I0,'C',C);
end
