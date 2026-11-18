function S=viii_storage_matrices(C,P8)
%VIII_STORAGE_MATRICES Local quadratic subsystem storage matrices.
%
% V_i = 0.5*y_i'*P_i*y_i in normalized reduced coordinates.
% F and R use the Hessians of the exact SO(3) storages at the origin.

B=viii_blocks();
P=cell(8,1);

P{1}=eye(3);
P{2}=eye(3);
P{3}=eye(3);
P{4}=eye(3);

% Governor: normalized second-order nominal filter Lyapunov matrix.
AG=[zeros(3),eye(3);
   -(C.R.Omega^2),-2*C.R.Zeta*C.R.Omega];
dG=P8.scale.G;
AGb=diag(1./dG)*AG*diag(dG);
P{5}=rg_solve_lyap(AGb,eye(6));

% Command-filter exact-storage Hessian.
HF=[C.F.kRf*eye(3), C.F.cF*eye(3);
    C.F.cF*eye(3), eye(3)];
dF=P8.scale.F;
P{6}=diag(dF)*HF*diag(dF);

% Rotational exact-storage Hessian.
I=C.P.I;
HR=[C.Dp.kR*eye(3), C.Dp.cR*I;
    C.Dp.cR*I, I];
dR=P8.scale.R;
P{7}=diag(dR)*HR*diag(dR);

% Actuator storage in normalized error coordinates.
P{8}=eye(4);

m=zeros(8,1); M=zeros(8,1);
for i=1:8
    ev=eig(0.5*(P{i}+P{i}.'));
    m(i)=0.5*min(ev);
    M(i)=0.5*max(ev);
end

S.P=P;
S.m=m;
S.M=M;
S.blocks=B;
end
