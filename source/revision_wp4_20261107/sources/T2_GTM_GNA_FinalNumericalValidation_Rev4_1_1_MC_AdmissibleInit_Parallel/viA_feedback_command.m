function C = viA_feedback_command(y,xf,R,A)
%VIA_FEEDBACK_COMMAND Section VI-A feedback force/path command.
%
% y  = [V;gamma;chi] actual air-path state
% xf = [Vf;gammaf;chif;nuV;nuGamma;nuChi] governed reference
%
% Ideal VI-A realization assumes beta=0 and alpha=alpha_c.

V=y(1); gamma=y(2); chi=y(3);
Vf=xf(1); gammaf=xf(2); chif=xf(3);
nuV=xf(4); nuG=xf(5); nuC=xf(6);

eV=V-Vf;
eGamma=gamma-gammaf;
eChi=rg_angle_error(chi,chif);

Xc=R.P.m*(nuV-A.kV*eV+R.g*sin(gamma));
aGamma=V*(nuG-A.kGamma*eGamma)+R.g*cos(gamma);
aChi=V*cos(gamma)*(nuC-A.kChi*eChi);
Nc=R.P.m*hypot(aGamma,aChi);
mu=atan2(aChi,aGamma);

[alpha,ri]=rg_alpha_root(V,Xc,Nc,R);
if ~ri.bracket
    error('viA_feedback_command:root','Feedback AoA root is not bracketed.');
end

Aa=rg_nominal_aero(V,alpha,R);
Fparallel=Aa.D+Xc;

% VI-A is the ideal coordinated case: beta=0 and alpha=alpha_c.
beta=0;
Tstar=Fparallel/(cos(alpha)*cos(beta));

C.e=[eV;eGamma;eChi];
C.X=Xc;
C.aGamma=aGamma;
C.aChi=aChi;
C.N=Nc;
C.mu=mu;
C.alpha=alpha;
C.beta=beta;
C.Fparallel=Fparallel;
C.Tstar=Tstar;
C.D=Aa.D;
C.L=Aa.L;
C.rootInfo=ri;
C.kappaT=cos(alpha)*cos(beta);
end
