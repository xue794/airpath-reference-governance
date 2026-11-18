function R = simulate_smc_guidance_case(mode,T,G,E)
%SIMULATE_SMC_GUIDANCE_CASE Integrate one isolated Section IV guidance case.

[rd0,~] = smc_guidance_reference(0);
r0 = rd0 + G.e0;
z0 = [r0; G.ze0; E.what0];

tg = linspace(0,T,round(20*T)+1);
opts = odeset('RelTol',G.RelTol,'AbsTol',G.AbsTol);
[t,z] = ode45(@(tt,zz)smc_guidance_augmented_rhs(tt,zz,G,E,mode), ...
              tg,z0,opts);

n = numel(t);
er=zeros(n,3); s=zeros(n,3); what=zeros(n,3); tw=zeros(n,3);
eg=zeros(n,3); ev=zeros(n,3); vc=zeros(n,3); vf=zeros(n,3); va=zeros(n,3);
w=zeros(n,3); nw=zeros(n,3); rd=zeros(n,3);

for k=1:n
    Y=smc_guidance_outputs(t(k),z(k,:).',G,E,mode);
    er(k,:)=Y.er.';
    s(k,:)=Y.s.';
    what(k,:)=z(k,7:9);
    tw(k,:)=Y.tildew.';
    eg(k,:)=Y.eg.';
    ev(k,:)=Y.ev.';
    vc(k,:)=Y.vc.';
    vf(k,:)=Y.vf.';
    va(k,:)=Y.va.';
    w(k,:)=Y.w.';
    nw(k,:)=Y.nw.';
    rd(k,:)=Y.rd.';
end

R=struct('t',t,'z',z,'r',z(:,1:3),'rd',rd,'er',er,'s',s, ...
         'what',what,'tildew',tw,'eg',eg,'ev',ev,'vc',vc,'vf',vf, ...
         'va',va,'w',w,'nw',nw,'mode',mode);
end
