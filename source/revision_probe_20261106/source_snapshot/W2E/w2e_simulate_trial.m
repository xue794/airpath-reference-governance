function S=w2e_simulate_trial(method,C,O)
%W2E_SIMULATE_TRIAL One full 120-s randomized trajectory using a pre-applied sample.
mode=O.Mode;K=cl6_case_info(mode,C);[x0,I0]=cl6_initial_state(mode,C);
xi0=rg_backup_state(x0(20:25),C.R);VB0=xi0.'*C.R.PB*xi0;
I0.VB0=VB0;I0.VB0_over_rho=VB0/C.R.rhoB;
tgrid=(0:C.sample_dt:K.T).';if tgrid(end)<K.T,tgrid=[tgrid;K.T];end
opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep);
wall=tic;[t,x]=ode45(@(tt,xx)w2e_rhs(tt,xx,mode,C,method,O),tgrid,x0,opts);integrationWall=toc(wall);
if any(~isfinite(x(:))),error('W2E:Nonfinite','Nonfinite state encountered.');end
outWall=tic;Y=w2e_outputs(t,x,mode,C,method,O);outputWall=toc(outWall);
S=struct('t',t,'x',x,'Y',Y,'mode',mode,'Method',method,'Case',K,'I0',I0,'C',C, ...
 'IntegrationWall_s',integrationWall,'OutputWall_s',outputWall);
end
