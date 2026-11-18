function S = simulate_viA_case(mode,T,R,A,E)
%SIMULATE_VIA_CASE Integrate Section V governor + VI-A ideal path dynamics.
z0=viA_initial_state(mode,R,A);
tg=linspace(0,T,round(T/A.sample_dt)+1);
opts=odeset('RelTol',A.RelTol,'AbsTol',A.AbsTol);
[t,z]=ode45(@(tt,zz)viA_combined_rhs(tt,zz,mode,R,A,E),tg,z0,opts);

n=numel(t);
e=zeros(n,3); alpha=zeros(n,1); Tstar=zeros(n,1); mu=zeros(n,1);
Fparallel=zeros(n,1); kappa=zeros(n,1);
aff=zeros(n,1); Tff=zeros(n,1); muff=zeros(n,1);
lambda=zeros(n,1); w=zeros(n,3); wdot=zeros(n,3);
cc=zeros(n,3);

for k=1:n
    xf=z(k,1:6).'; y=z(k,7:9).';
    C=viA_feedback_command(y,xf,R,A);
    F=rg_feedforward(xf,R);
    sig=viA_case_signal(t(k),mode,E);

    [vc,~]=rg_raw_test_command(t(k),sig.refMode,R);
    [cck,~]=rg_raw_command_chart(vc,R);
    u0=rg_nominal_control(xf,cck,R);
    O=rg_scalar_governor(xf,u0,R);

    e(k,:)=C.e.';
    alpha(k)=C.alpha;
    Tstar(k)=C.Tstar;
    mu(k)=C.mu;
    Fparallel(k)=C.Fparallel;
    kappa(k)=C.kappaT;

    aff(k)=F.alpha;
    Tff(k)=F.T;
    muff(k)=F.mu;
    lambda(k)=O.lambda;
    w(k,:)=sig.w.';
    wdot(k,:)=sig.wdot.';
    cc(k,:)=cck.';
end

S=struct('t',t,'z',z,'xf',z(:,1:6),'y',z(:,7:9),'e',e, ...
         'alpha',alpha,'Tstar',Tstar,'mu',mu,'Fparallel',Fparallel, ...
         'kappa',kappa,'alphaFF',aff,'Tff',Tff,'muFF',muff, ...
         'lambda',lambda,'w',w,'wdot',wdot,'cc',cc,'mode',mode);
end
