function preflight_comparison()
here=fileparts(mfilename('fullpath'));src=fullfile(fileparts(here),'revision_probe_20261106','source_snapshot');
addpath(fullfile(src,'core'));addpath(fullfile(src,'W2D'));addpath(fullfile(here,'full6dof'));addpath(fullfile(here,'aoa'));
R=rg_params();R.rhoB=.1;
states={[31;0;0;0;0;0],[38.1;.02;.1;.2;.003;.01]};
controls={[.2;.02;.01],[-.4;-.01;.015]};maxError=0;
for i=1:numel(states)
 x=states{i};u=controls{i};v=[x(4:6);u];C=cp6_affine(x,R);
 for mode={'A0','Am'}
  A=cp_aoa_rows(x,R,mode{1},C.alpha);dt=1e-5;
  Ap=cp_aoa_rows(x+dt*v,R,mode{1});Am=cp_aoa_rows(x-dt*v,R,mode{1});
  analytic=A.c+A.G*u-R.kN*A.h;numeric=(Ap.h-Am.h)/(2*dt);
  maxError=max(maxError,max(abs(analytic-numeric)));
  assert(max(abs(analytic-numeric))<1e-7,'CP:Derivative','Angle margin derivative mismatch');
 end
end
fprintf('CP_PREFLIGHT_PASS max_derivative_error=%.12g rad/s\n',maxError);
end
