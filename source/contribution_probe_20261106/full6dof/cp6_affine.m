function C = cp6_affine(x,R)
%RG_CONSTRAINT_AFFINE Construct phi_j(x,u)=c_j(x)+g_j(x)'u >= 0.
% Includes HOCBF state boxes, rate CBFs, normal-force, bank, thrust,
% positive a_gamma, and backup-ellipsoid CBF.

V=x(1); gam=x(2); nuV=x(4); nuG=x(5); nuC=x(6);
m=R.P.m; gg=R.g;
names={}; cc=[]; GG=[];

    function add(name,c,grow)
        names{end+1,1}=name; %#ok<AGROW>
        cc(end+1,1)=c; %#ok<AGROW>
        GG(end+1,:)=grow(:).'; %#ok<AGROW>
    end

% Relative-degree-two speed and climb constraints.
add('V_min',(R.k1+R.k2)*nuV+R.k1*R.k2*(V-R.Vmin),[1 0 0]);
add('V_max',-(R.k1+R.k2)*nuV+R.k1*R.k2*(R.Vmax-V),[-1 0 0]);
add('gamma_min',(R.k1+R.k2)*nuG+R.k1*R.k2*(gam+R.gammaMax),[0 1 0]);
add('gamma_max',-(R.k1+R.k2)*nuG+R.k1*R.k2*(R.gammaMax-gam),[0 -1 0]);

% First-order rate CBFs.
nu=[nuV;nuG;nuC];
for i=1:3
    ei=zeros(1,3); ei(i)=1;
    add(sprintf('nu%d_min',i),R.kRate*(nu(i)+R.nuMax(i)),ei);
    add(sprintf('nu%d_max',i),R.kRate*(R.nuMax(i)-nu(i)),-ei);
end

% Force geometry and its affine derivatives.
ag=V*nuG+gg*cos(gam);
ac=V*cos(gam)*nuC;
dag_c=nuV*nuG-gg*sin(gam)*nuG;
dag_g=[0 V 0];
dac_c=nuV*cos(gam)*nuC-V*sin(gam)*nuG*nuC;
dac_g=[0 0 V*cos(gam)];

add('aGamma_min',dag_c+R.kAGamma*(ag-R.aGammaMin),dag_g);

tm=tan(R.muMax);
hmp=tm*ag-ac;
hmm=tm*ag+ac;
add('mu_plus',tm*dag_c-dac_c+R.kMu*hmp,tm*dag_g-dac_g);
add('mu_minus',tm*dag_c+dac_c+R.kMu*hmm,tm*dag_g+dac_g);

X=m*(nuV+gg*sin(gam));
Xdot_c=m*gg*cos(gam)*nuG;
Xdot_g=[m 0 0];
rr=hypot(ag,ac);
N=m*rr;
Ndot_c=m*(ag*dag_c+ac*dac_c)/rr;
Ndot_g=m*(ag*dag_g+ac*dac_g)/rr;

% Normal-force barriers at the certified AoA endpoints.
Mlo=rg_N_map(R.alphaLo,V,X,R);
Mhi=rg_N_map(R.alphaHi,V,X,R);
MloDot_c=Mlo.V*nuV+Mlo.X*Xdot_c;
MloDot_g=Mlo.X*Xdot_g;
MhiDot_c=Mhi.V*nuV+Mhi.X*Xdot_c;
MhiDot_g=Mhi.X*Xdot_g;

hNm=N-Mlo.value-R.DeltaN;
hNp=Mhi.value-N-R.DeltaN;
add('N_minus',Ndot_c-MloDot_c+R.kN*hNm,Ndot_g-MloDot_g);
add('N_plus',MhiDot_c-Ndot_c+R.kN*hNp,MhiDot_g-Ndot_g);

% Implicit AoA and thrust derivative, separated into constant and u terms.
[alpha,ri]=cp6_alpha_root(V,X,N,R,'reference');
if ~ri.bracket
    error('rg_constraint_affine:root','AoA inversion bracket lost.');
end
Ma=rg_N_map(alpha,V,X,R);

adot_c=(Ndot_c-Ma.V*nuV-Ma.X*Xdot_c)/Ma.alpha;
adot_g=(Ndot_g-Ma.X*Xdot_g)/Ma.alpha;

A=rg_nominal_aero(V,alpha,R);
Ddot_c=A.DV*nuV+A.Da*adot_c;
Ddot_g=A.Da*adot_g;
T=(A.D+X)/cos(alpha);
Tdot_c=(Ddot_c+Xdot_c)/cos(alpha)+(A.D+X)/cos(alpha)*tan(alpha)*adot_c;
Tdot_g=(Ddot_g+Xdot_g)/cos(alpha)+(A.D+X)/cos(alpha)*tan(alpha)*adot_g;

add('T_min',Tdot_c+R.kT*(T-(R.Tmin+R.DeltaT)),Tdot_g);
add('T_max',-Tdot_c+R.kT*((R.Tmax-R.DeltaT)-T),-Tdot_g);

% Backup ellipsoid CBF.
xi=rg_backup_state(x,R);
y=R.PB*xi;
hB=R.rhoB-xi'*y;
hBdot_c=-2*(y(1)*nuV+y(2)*nuG);
hBdot_g=-2*[y(3) y(4) y(5)];
add('backup',hBdot_c+R.kB*hB,hBdot_g);

C.c=cc;
C.G=GG;
C.names=names;
C.alpha=alpha;
C.T=T;
C.N=N;
C.aGamma=ag;
C.aChi=ac;
C.mu=atan2(ac,ag);
C.hNminus=hNm;
C.hNplus=hNp;
C.hB=hB;

% HOCBF initialization variables.
C.h=[V-R.Vmin;
     R.Vmax-V;
     gam+R.gammaMax;
     R.gammaMax-gam];
C.psi1=[nuV+R.k1*(V-R.Vmin);
       -nuV+R.k1*(R.Vmax-V);
        nuG+R.k1*(gam+R.gammaMax);
       -nuG+R.k1*(R.gammaMax-gam)];
end
