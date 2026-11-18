function S=w2d_simulate_case(mode,method,C,O)
%W2D_SIMULATE_CASE One targeted full-6DOF trajectory from the certified trim.
K=w2d_case_info(mode,C,O);
% The target profiles exactly equal the source hold reference at t=0 and use still wind.
% Therefore the frozen hold initializer provides the same certified trim without modifying source code.
[x0,I0]=cl6_initial_state('hold',C);[rd0,rdDot0]=w2d_reference(0,mode,C,O);
assert(norm(I0.rd-rd0,inf)<1e-12 && norm(I0.rdDot-rdDot0,inf)<1e-12,'W2D:InitialReference','Target profile does not match frozen trim at t=0.');
xi0=rg_backup_state(x0(20:25),C.R);VB0=xi0.'*C.R.PB*xi0;assert(VB0<=C.R.rhoB+O.InitialBackupTol,'W2D:InitialCertificate','Target case starts outside K_B.');
tgrid=(0:C.sample_dt:K.T).';if tgrid(end)<K.T,tgrid=[tgrid;K.T];end
opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep);wall=tic;[t,x]=ode45(@(tt,xx)w2d_rhs(tt,xx,mode,C,method,O),tgrid,x0,opts);integrationWall=toc(wall);
outWall=tic;Y=w2d_outputs(t,x,mode,C,method,O);outputWall=toc(outWall);
I0.TargetMode=mode;I0.VB0=VB0;I0.VB0_over_rho=VB0/C.R.rhoB;I0.TargetAirPath=K.TargetAirPath;
S=struct('t',t,'x',x,'Y',Y,'mode',mode,'Method',method,'Case',K,'I0',I0,'C',C,'IntegrationWall_s',integrationWall,'OutputWall_s',outputWall);
end
