function audit_alpha_interpretation()
% Same-state interpretation only; no new closed-loop comparison or tuning.
here=fileparts(mfilename('fullpath'));
base='D:/Research/paper1_Submission/revision_probe_20261106/source_snapshot';
addpath([base '/core']);addpath([base '/W2E']);addpath(here);
R=rg_params();R.rhoB=.10;Q=w2b_options();Q.FeasibleStart=zeros(3,1);
x=[31;0;0;0;0;0];request=[31;8*pi/180;0];u0=rg_nominal_control(x,request,R);
C10=rg_constraint_affine(x,R);R0=R;R0.DeltaN=0;C0=rg_constraint_affine(x,R0);
F=rg_feedforward(x,R);Ma=rg_N_map(F.alpha,x(1),F.X,R);
Mhi=rg_N_map(R.alphaHi,x(1),F.X,R);
names={'thrust_only','endpoint_DeltaN0','endpoint_DeltaN10'};
sets={[1:13 16:17],1:17,1:17};allC={C10,C0,C10};Rows=struct([]);
for j=1:3
 C=allC{j};ix=sets{j};[u,D]=w2b_qp3(u0,-C.G(ix,:),C.c(ix),Q.Scale,Q);assert(D.Success);
 % Same analytical chain rule as the frozen kernel, kept explicit for audit.
 V=x(1);gam=x(2);nuV=x(4);nuG=x(5);nuChi=x(6);m=R.P.m;
 ag=V*nuG+R.g*cos(gam);ac=V*cos(gam)*nuChi;
 dag=nuV*nuG-R.g*sin(gam)*nuG+V*u(2);
 dac=nuV*cos(gam)*nuChi-V*sin(gam)*nuG*nuChi+V*cos(gam)*u(3);
 Ndot=m*(ag*dag+ac*dac)/hypot(ag,ac);
 Xdot=m*(R.g*cos(gam)*nuG+u(1));
 alphadot=(Ndot-Ma.V*nuV-Ma.X*Xdot)/Ma.alpha;
 gap=R.alphaHi-F.alpha;
 anglePhi=-alphadot+R.kN*gap;
 endpointPhi0=C0.c(15)+C0.G(15,:)*u;
 endpointPhi10=C10.c(15)+C10.G(15,:)*u;
 hNraw=F.hNplus+R.DeltaN;
 hNdot=Mhi.V*nuV+Mhi.X*Xdot-Ndot;
 eta=hNraw/gap;
 etaDot=(hNdot*gap+hNraw*alphadot)/gap^2;
 identity=eta*(-alphadot+(R.kN+etaDot/eta)*gap);
 q=struct('method',names{j},'u',u','alpha_dot_radps',alphadot,'alpha_dot_degps',alphadot*180/pi, ...
  'standard_alpha_upper_CBF_radps',anglePhi,'standard_alpha_upper_CBF_degps',anglePhi*180/pi, ...
  'endpoint_DeltaN0_CBF_Nps',endpointPhi0,'endpoint_DeltaN10_CBF_Nps',endpointPhi10, ...
  'eta_N_per_rad',eta,'eta_dot_N_per_rad_s',etaDot, ...
  'state_dependent_equivalent_angle_gain_per_s',R.kN+etaDot/eta, ...
  'identity_error_Nps',identity-endpointPhi0, ...
  'thrust_lower_CBF_Nps',C10.c(16)+C10.G(16,:)*u, ...
  'thrust_upper_CBF_Nps',C10.c(17)+C10.G(17,:)*u, ...
  'maxKKT',max([D.PrimalNormalized D.StationarityNormalized D.ComplementarityNormalized]));
 Rows=[Rows;q]; %#ok<AGROW>
end
out=struct('state',x','request',request','alpha_deg',F.alpha*180/pi, ...
 'alpha_upper_deg',R.alphaHi*180/pi,'angular_gap_deg',(R.alphaHi-F.alpha)*180/pi, ...
 'Tff_N',F.T,'raw_endpoint_margin_N',F.hNplus+R.DeltaN, ...
 'standard_alpha_gain_per_s',R.kN,'rows',Rows, ...
 'scope','Fixed-state analytic interpretation; no alpha-CBF baseline trajectory was run');
fid=fopen(fullfile(here,'alpha_interpretation.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);
disp(out);disp(Rows);fprintf('ALPHA_INTERPRETATION_DONE\n');
end
