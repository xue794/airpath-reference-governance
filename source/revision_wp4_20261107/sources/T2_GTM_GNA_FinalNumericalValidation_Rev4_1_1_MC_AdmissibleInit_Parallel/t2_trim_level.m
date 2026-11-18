function TR = t2_trim_level(V,h,P)
%T2_TRIM_LEVEL  Straight/level trim for the reduced point-thrust GNA plant.
% Unknowns: alpha, elevator, thrust. beta=phi=gamma=p=q=r=da=dr=0.
% Uses base-MATLAB fminsearch; no Optimization Toolbox required.

rho=t2_isa_density(h);
y0=[deg2rad(5);deg2rad(1);18.5];
scale=[P.g;P.g;P.g*P.cbar];
obj=@(y) sum((residual(y)./scale).^2);
opts=optimset('TolX',1e-12,'TolFun',1e-14,'MaxIter',3000,'MaxFunEvals',10000,'Display','off');
y=fminsearch(obj,y0,opts);
R=residual(y);
TR=struct('alpha',y(1),'deltae',y(2),'T',y(3),'V',V,'h',h, ...
          'residual',R,'objective',obj(y));

    function R=residual(y)
        alpha=y(1); de=y(2); T=y(3);
        vgB=V*[cos(alpha);0;sin(alpha)];
        CIB=euler321_dcm(0,alpha,0);  % gamma=0 => theta=alpha
        [fA,tauA]=t2_gna_aero(vgB,[0;0;0],[0;de;0],CIB,[0;0;0],h,P);
        zeta=CIB.'*[0;0;1];
        aB=P.g*zeta+(fA+[T;0;0])/P.m;
        alpha_ddot_equiv = tauA(2)/(P.m*P.cbar); % only for residual scaling
        R=[aB(1);aB(3);alpha_ddot_equiv];
    end
end
