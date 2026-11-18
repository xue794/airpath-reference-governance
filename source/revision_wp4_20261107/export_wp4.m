function export_wp4()
% Reconstruct saved samples with the frozen MATLAB controller; never integrate.
root=fileparts(mfilename('fullpath')); old=path; cleanup=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(root,'sources','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel'),'-begin');
addpath(fullfile(root,'sources','P1_Aerospace_W2C_R2018a_v1_0'),'-begin');
addpath(fullfile(root,'sources','transparent'),'-begin');
rows=jsondecode(fileread(fullfile(root,'inputs','new_trajectories.json'))); C0=cl6_params();C0.R.rhoB=.10;
out=fullfile(root,'results'); O=w2c_options();
fid=fopen(fullfile(root,'inputs','configuration.json'),'w');fprintf(fid,'%s',jsonencode(C0));fclose(fid);
Summaries={};
for j=1:numel(rows)
 r=rows(j);L=load(r.path);if isfield(L,'S'),B=L.S;else,B=L.Baseline;end
 storedRho=B.C.R.rhoB;C=B.C;C.R.rhoB=.10;assert(isequaln(C,C0),'WP3:Config','Configuration differs.');
 assert(strcmp(B.mode,r.case));t=B.t(:);x=B.x;n=numel(t);
 assert(size(x,2)==36 && n==size(x,1) && t(end)<=r.duration_s+1e-10);
 assert(t(1)==0 && max(abs(diff(t)-.05))<1e-10);
 [x0,~]=cl6_initial_state(B.mode,C);initError=max(abs(x(1,:).'-x0));assert(initError<1e-12);
 names={'time_s','Vf','gammaf','chif','nuV','nuGamma','nuChi','eV','eGamma','eChi','thetaC','thetaF','thetaR', ...
 'ref_ag','ref_ac','ref_X','ref_N','ref_lower','ref_upper','ref_T','ref_alpha','ref_hB', ...
 'fb_ag','fb_ac','fb_X','fb_N','fb_lower','fb_upper','Tstar','Traw','Tact','kappa','alphaC','muC', ...
 'actual_V','actual_gamma','actual_chi','actual_alpha','actual_beta','position_error', ...
 'wind_N','wind_E','wind_D','correction','intervened','lambda'};
 D=zeros(n,numel(names));Ref=zeros(n,7);Feed=zeros(n,7);Phys=zeros(n,16);Domain=zeros(n,6);Phi=zeros(n,18);
 identity=zeros(n,5);projected=zeros(n,1);savedFields={'V','Tcmd','thetaC'};savedAvailable=cellfun(@(f)isfield(B.Y,f),savedFields);
 for k=1:n
  xx=x(k,:).';xf=xx(20:25);R=C.R;
  if strcmp(r.method,'Transparent'),S=cl6_control_snapshot_transparent_R2018a(t(k),xx,B.mode,C);
  else,S=w2c_control_snapshot(t(k),xx,B.mode,C,r.method,O);end
  A=S.Gov.C;F=S.Cmd;Air=S.Air;
  ml=rg_N_map(R.alphaLo,Air.V,F.X,R);mh=rg_N_map(R.alphaHi,Air.V,F.X,R);
  b=[F.N-ml.value;mh.value-F.N];rfb=[A.hNminus;A.hNplus]+R.DeltaN;
  Xf=R.P.m*(xf(4)+R.g*sin(xf(2)));af=rg_N_map(A.alpha,xf(1),Xf,R);
  thetaC=dcm_distance(S.Cc,Air.CIB);thetaF=dcm_distance(S.Cc,S.Cr);thetaR=dcm_distance(S.Cr,Air.CIB);
  if strcmp(r.method,'Transparent'),correction=0;intervened=0;lam=1;
  else,correction=S.Gov.CorrectionNorm;intervened=S.Gov.Intervened;if isfield(S.Gov,'lambda'),lam=S.Gov.lambda;else,lam=NaN;end,end
  D(k,:)=[t(k);xf;F.e;thetaC;thetaF;thetaR;A.aGamma;A.aChi;Xf;A.N;rfb;A.T;A.alpha;A.hB; ...
   F.aGamma;F.aChi;F.X;F.N;b;F.Tstar;S.Tcmd;xx(36);Air.kappa;F.alpha;F.mu; ...
   Air.V;Air.gamma;Air.chi;Air.alpha;Air.beta;norm(S.Gd.er);Air.wI;correction;intervened;lam].';
  Ref(k,:)=[A.aGamma-R.aGammaMin;tan(R.muMax)*A.aGamma-A.aChi;tan(R.muMax)*A.aGamma+A.aChi; ...
    rfb;A.T-R.Tmin;R.Tmax-A.T].';
  Feed(k,:)=[F.aGamma-R.aGammaMin;tan(R.muMax)*F.aGamma-F.aChi;tan(R.muMax)*F.aGamma+F.aChi;b;S.Tcmd-R.Tmin;R.Tmax-S.Tcmd].';
  P=C.P7;Phys(k,:)=[xx(33:35)-P.deltaMin;P.deltaMax-xx(33:35);xx(36)-P.Tmin;P.Tmax-xx(36); ...
    S.deltaDot+P.deltaRateMax;P.deltaRateMax-S.deltaDot;S.TDot+P.TRateMax;P.TRateMax-S.TDot].';
  Domain(k,:)=[Air.V-R.Vmin;R.Vmax-Air.V;Air.gamma+R.gammaMax;R.gammaMax-Air.gamma;Air.alpha-R.alphaPhysicalMin;R.alphaPhysicalMax-Air.alpha].';
  Phi(k,:)=S.Gov.phi.';projected(k)=S.Tcmd<P.Tmin || S.Tcmd>P.Tmax;
  identity(k,1:2)=[abs(A.N-af.value),max(abs(S.Gd.er.'-B.Y.er(k,:)))];
  checkValues=[Air.V,S.Tcmd,thetaC];
  for ii=1:3,if savedAvailable(ii),identity(k,ii+2)=abs(checkValues(ii)-B.Y.(savedFields{ii})(k));end,end
 end
 assert(max(max(identity(:,2:5)))<1e-9);
 key=[r.method '_' r.case];
 save(fullfile(out,[key '_export.mat']),'t','x','D','names','C','Ref','Feed','Phys','Domain','Phi','-v7');
 writetable(array2table(D,'VariableNames',names),fullfile(out,[key '_samples.csv']));
 Z=struct('Case',r.case,'Method',r.method,'Samples',n,'End_s',t(end),'InitialState',x0,'InitError',initError, ...
  'StoredRho',storedRho,'AuditRho',C.R.rhoB,'ConfigurationMatch',true,'MaxSavedIdentityErrors',max(identity,[],1), ...
  'OptionalSavedFields',{savedFields},'OptionalSavedFieldAvailable',savedAvailable, ...
  'PositionRMS',sqrt(mean(D(:,40).^2)),'PositionPeak',max(D(:,40)), ...
  'ReferenceMinimumMargins',min(Ref,[],1),'FeedbackMinimumMargins',min(Feed,[],1), ...
  'ReferenceViolationCounts',sum(Ref < -1e-8,1),'FeedbackViolationCounts',sum(Feed < -1e-8,1), ...
  'PhysicalActuatorViolationSamples',sum(any(Phys < -1e-8,2)),'StudyDomainViolationSamples',sum(any(Domain < -1e-8,2)), ...
  'ThrustProjectionSamples',sum(projected),'Source',r.path);
 Summaries{j}=Z; %#ok<AGROW>
 fprintf('EXPORTED %s %d samples; saved identity %.3g\n',key,n,max(max(identity(:,2:5))));
end
Summary=[Summaries{:}];fid=fopen(fullfile(out,'export_summary.json'),'w');fprintf(fid,'%s',jsonencode(Summary));fclose(fid);
fprintf('WP4_EXPORT_COMPLETE no integrations\n');
end
