function S = simulate_reference_governor(mode,T,R)
%SIMULATE_REFERENCE_GOVERNOR Integrate Section V governor only.
if strcmpi(mode,'wrap')
    [~,c0]=rg_raw_test_command(0,'wrap',R);
    chi0=c0(3);
else
    chi0=0;
end
x0=[R.Vb;R.gammaB;chi0;0;0;0];
tg=linspace(0,T,round(T/R.sample_dt)+1);
opts=odeset('RelTol',R.RelTol,'AbsTol',R.AbsTol);
[t,x]=ode45(@(tt,xx)rg_governor_rhs(tt,xx,mode,R),tg,x0,opts);

n=numel(t);
cc=zeros(n,3); u=zeros(n,3); u0=zeros(n,3); lam=zeros(n,1);
minphi=zeros(n,1); alpha=zeros(n,1); Tff=zeros(n,1); mu=zeros(n,1);
hB=zeros(n,1); hN=zeros(n,2); etaChi=zeros(n,1);

for k=1:n
    [vc,~]=rg_raw_test_command(t(k),mode,R);
    [cck,~]=rg_raw_command_chart(vc,R);
    u0k=rg_nominal_control(x(k,:).',cck,R);
    O=rg_scalar_governor(x(k,:).',u0k,R);
    cc(k,:)=cck.';
    u(k,:)=O.u.';
    u0(k,:)=u0k.';
    lam(k)=O.lambda;
    minphi(k)=min(O.phi);
    alpha(k)=O.C.alpha;
    Tff(k)=O.C.T;
    mu(k)=O.C.mu;
    hB(k)=O.C.hB;
    hN(k,:)=[O.C.hNminus O.C.hNplus];
    etaChi(k)=rg_angle_error(x(k,3),cck(3));
end

S=struct('t',t,'x',x,'cc',cc,'u',u,'u0',u0,'lambda',lam, ...
         'minphi',minphi,'alpha',alpha,'Tff',Tff,'mu',mu,'hB',hB, ...
         'hN',hN,'etaChi',etaChi,'mode',mode);
end
