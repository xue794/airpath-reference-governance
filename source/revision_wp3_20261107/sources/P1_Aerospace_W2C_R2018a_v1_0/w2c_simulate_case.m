function S=w2c_simulate_case(mode,method,C,O)
%W2C_SIMULATE_CASE Integrate one deterministic full-6DOF comparison trajectory.
K=cl6_case_info(mode,C);[x0,I0]=cl6_initial_state(mode,C);
tgrid=(0:C.sample_dt:K.T).';if tgrid(end)<K.T,tgrid=[tgrid;K.T];end
opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep);
wall=tic;[t,x]=ode45(@(tt,xx)w2c_rhs(tt,xx,mode,C,method,O),tgrid,x0,opts);integrationWall=toc(wall);
outWall=tic;Y=w2c_outputs(t,x,mode,C,method,O);outputWall=toc(outWall);
S=struct('t',t,'x',x,'Y',Y,'mode',mode,'Method',method,'Case',K,'I0',I0,'C',C, ...
    'IntegrationWall_s',integrationWall,'OutputWall_s',outputWall);
end
