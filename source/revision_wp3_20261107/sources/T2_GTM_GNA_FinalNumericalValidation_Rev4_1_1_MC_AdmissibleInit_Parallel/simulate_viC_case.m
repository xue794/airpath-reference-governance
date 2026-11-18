function S=simulate_viC_case(mode,T,R,A,F,Ew)
%SIMULATE_VIC_CASE Section V/VI-A upstream + dynamic SO(3) command filter.
%
% The upstream reference/path dynamics are independent of the command filter
% in VI-C.  They are therefore generated once, then used as a smooth pchip
% command input to the filter.  omega_c is evaluated analytically at the
% stored upstream samples for diagnostics and error variables.

Up=simulate_viA_case(mode,T,R,A,Ew);

Cmd.t=Up.t;
Cmd.gamma=Up.y(:,2);
Cmd.chi=Up.y(:,3);
Cmd.mu=Up.mu;
Cmd.alpha=Up.alpha;

Cc0=viB_command_attitude(Cmd.gamma(1),Cmd.chi(1),Cmd.mu(1),Cmd.alpha(1));
qr0=viC_dcm_to_quat(Cc0);
x0=[qr0;zeros(3,1)];

opts=odeset('RelTol',F.RelTol,'AbsTol',F.AbsTol,'MaxStep',F.filterMaxStep);
[t,x]=ode45(@(tt,xx)viC_filter_interp_rhs(tt,xx,Cmd,F),Up.t,x0,opts);

n=numel(t);
theta=zeros(n,1);
Psi=zeros(n,1);
e=zeros(n,3);
xi=zeros(n,3);
omegaC=zeros(n,3);
omegaR=x(:,5:7);
qNorm=zeros(n,1);
Vstorage=zeros(n,1);
alphaDot=zeros(n,1);
muDot=zeros(n,1);

for k=1:n
    Cc=viC_interp_command(t(k),Cmd);
    Cr=quat_to_dcm(x(k,1:4).');

    % t uses exactly the same sample grid as Up.t.
    Kin=viC_command_kinematics(t(k),Up.z(k,:).',mode,R,A,Ew);
    Err=viC_filter_error(Cc,Kin.omegaC,Cr,omegaR(k,:).');

    theta(k)=Err.theta;
    Psi(k)=Err.Psi;
    e(k,:)=Err.e.';
    xi(k,:)=Err.xi.';
    omegaC(k,:)=Kin.omegaC.';
    qNorm(k)=norm(x(k,1:4));
    Vstorage(k)=0.5*(omegaR(k,:)*omegaR(k,:).') + ...
                F.kRf*Err.Psi + F.cF*(Err.e.'*omegaR(k,:).');
    alphaDot(k)=Kin.alphaDot;
    muDot(k)=Kin.muDot;
end

S=struct('t',t,'x',x,'Up',Up,'Cmd',Cmd,'theta',theta,'Psi',Psi, ...
         'e',e,'xi',xi,'omegaC',omegaC,'omegaR',omegaR, ...
         'qNorm',qNorm,'Vstorage',Vstorage,'alphaDot',alphaDot, ...
         'muDot',muDot,'mode',mode);
end
