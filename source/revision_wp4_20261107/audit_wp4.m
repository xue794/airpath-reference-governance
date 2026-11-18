function audit_wp4()
% Common layered audit on all nine original-case/method cells. No integrations.
H=fileparts(mfilename('fullpath'));addpath(fullfile(H,'sources','T2_GTM_GNA_FinalNumericalValidation_Rev4_1_1_MC_AdmissibleInit_Parallel'),'-begin');
addpath(fullfile(H,'sources','P1_Aerospace_W2C_R2018a_v1_0'),'-begin');addpath(fullfile(H,'sources','transparent'),'-begin');
records=jsondecode(fileread(fullfile(H,'inputs','all_trajectories.json')));O=w2c_options();C0=cl6_params();C0.R.rhoB=.10;
out=fullfile(H,'results');allRows={};allSummary={};
refnames={'V_min','V_max','gamma_min','gamma_max','psi1_V_min','psi1_V_max','psi1_gamma_min','psi1_gamma_max', ...
 'nuV_min','nuGamma_min','nuChi_min','nuV_max','nuGamma_max','nuChi_max','aGamma_min','bank_plus','bank_minus', ...
 'raw_endpoint_lower','raw_endpoint_upper','reserved_endpoint_lower','reserved_endpoint_upper','reference_T_0N','reference_T_60N','reference_T_3N','reference_T_57N','backup'};
refunits={'m/s','m/s','rad','rad','m/s2','m/s2','rad/s','rad/s','m/s2','rad/s','rad/s','m/s2','rad/s','rad/s','m/s2','m/s2','m/s2','N','N','N','N','N','N','N','N','weighted_set'};
feednames={'aGamma_min','bank_plus','bank_minus','raw_endpoint_lower','raw_endpoint_upper','reserved_endpoint_lower','reserved_endpoint_upper', ...
 'ideal_T_0N','ideal_T_60N','ideal_T_3N','ideal_T_57N','raw_T_0N','raw_T_60N','raw_T_3N','raw_T_57N'};
feedunits=[repmat({'m/s2'},1,3),repmat({'N'},1,12)];
physnames={'aileron_min','elevator_min','rudder_min','aileron_max','elevator_max','rudder_max','actuator_T_min','actuator_T_max', ...
 'aileron_rate_min','elevator_rate_min','rudder_rate_min','aileron_rate_max','elevator_rate_max','rudder_rate_max','T_rate_min','T_rate_max'};
physunits=[repmat({'rad'},1,6),{'N','N'},repmat({'rad/s'},1,6),{'N/s','N/s'}];
domnames={'actual_V_min_study','actual_V_max_study','actual_gamma_min_study','actual_gamma_max_study','actual_alpha_min_study','actual_alpha_max_study'};
valueNames={'time_s','position_error_m','reference_T_N','ideal_T_N','raw_T_N','actuator_T_N','VB_over_rho', ...
 'reference_alpha_deg','command_alpha_deg','actual_alpha_deg','actual_beta_deg','kappa','thetaF_deg','thetaR_deg','thetaC_deg', ...
 'raw_chart_speed_margin','raw_chart_horizontal_margin','raw_phi_violation','normalized_phi_violation','governor_correction_norm','governor_intervened', ...
 'scalar_lambda','thrust_target_clipped','thrust_rate_limited','surface_rate_limited','allocator_bound_active','allocation_residual_norm','servo_residual_norm'};
for j=1:numel(records)
 if iscell(records),r=records{j};else,r=records(j);end
 L=load(r.path);if isfield(L,'S'),B=L.S;else,B=L.Baseline;end
 C=B.C;storedRho=C.R.rhoB;C.R.rhoB=.10;assert(isequaln(C,C0));t=B.t(:);x=B.x;n=numel(t);
 [x0,~]=cl6_initial_state(r.case,C);initDiff=max(abs(x0-x(1,:).'));assert(initDiff<1e-12);
 Ref=zeros(n,26);Feed=zeros(n,15);Phys=zeros(n,16);Domain=zeros(n,6);Phi=zeros(n,18);Values=zeros(n,numel(valueNames));er=zeros(n,3);
 QPKKT=nan(n,4);qpFallback=false(n,1);u0=zeros(n,3);uApplied=zeros(n,3);qpaths=cell(n,1);referenceRoot=true(n,1);feedbackRoot=true(n,1);
 for k=1:n
  xx=x(k,:).';xf=xx(20:25);R=C.R;P=C.P7;
  if strcmp(r.method,'Transparent'),S=cl6_control_snapshot_transparent_R2018a(t(k),xx,r.case,C);
  else,S=w2c_control_snapshot(t(k),xx,r.case,C,r.method,O);end
  A=S.Gov.C;F=S.Cmd;nu=xf(4:6);rb=[A.hNminus;A.hNplus]+R.DeltaN;
  lo=rg_N_map(R.alphaLo,S.Air.V,F.X,R);hi=rg_N_map(R.alphaHi,S.Air.V,F.X,R);fb=[F.N-lo.value;hi.value-F.N];
  Ref(k,:)=[A.h;A.psi1;nu+R.nuMax;R.nuMax-nu;A.aGamma-R.aGammaMin;tan(R.muMax)*A.aGamma-A.aChi;tan(R.muMax)*A.aGamma+A.aChi; ...
   rb;rb-R.DeltaN;A.T-R.Tmin;R.Tmax-A.T;A.T-R.Tmin-R.DeltaT;R.Tmax-R.DeltaT-A.T;A.hB].';
  Feed(k,:)=[F.aGamma-R.aGammaMin;tan(R.muMax)*F.aGamma-F.aChi;tan(R.muMax)*F.aGamma+F.aChi;fb;fb-R.DeltaN; ...
   F.Tstar-R.Tmin;R.Tmax-F.Tstar;F.Tstar-R.Tmin-R.DeltaT;R.Tmax-R.DeltaT-F.Tstar; ...
   S.Tcmd-R.Tmin;R.Tmax-S.Tcmd;S.Tcmd-R.Tmin-R.DeltaT;R.Tmax-R.DeltaT-S.Tcmd].';
  Phys(k,:)=[xx(33:35)-P.deltaMin;P.deltaMax-xx(33:35);xx(36)-P.Tmin;P.Tmax-xx(36); ...
   S.deltaDot+P.deltaRateMax;P.deltaRateMax-S.deltaDot;S.TDot+P.TRateMax;P.TRateMax-S.TDot].';
  Domain(k,:)=[S.Air.V-R.Vmin;R.Vmax-S.Air.V;S.Air.gamma+R.gammaMax;R.gammaMax-S.Air.gamma;S.Air.alpha-R.alphaPhysicalMin;R.alphaPhysicalMax-S.Air.alpha].';
  u0(k,:)=rg_nominal_control(xf,S.cc,R).';uApplied(k,:)=S.Gov.u.';Phi(k,:)=S.Gov.phi.';
  den=1+abs(A.c)+abs(A.G)*abs(S.Gov.u);normPhi=max([0;-S.Gov.phi./den]);
  if strcmp(r.method,'Transparent'),corr=0;intervened=0;lam=1;qpaths{k}='NA';
  else
   corr=S.Gov.CorrectionNorm;intervened=S.Gov.Intervened;
   if strcmp(r.method,'Proposed'),lam=S.Gov.lambda;qpaths{k}='NA';
   else,lam=NaN;Q=S.Gov.Solver;QPKKT(k,:)=[Q.PrimalNormalized,Q.StationarityNormalized,Q.DualViolationNormalized,Q.ComplementarityNormalized];qpFallback(k)=Q.FallbackUsed;qpaths{k}=Q.Path;end
  end
  targetT=min(max(S.Tcmd,P.Tmin),P.Tmax);preTdot=(targetT-xx(36))/P.tauT;
  preDelta=(S.deltaCmd-xx(33:35))./P.tauDelta;
  er(k,:)=S.Gd.er.';referenceRoot(k)=all(rb>=0);feedbackRoot(k)=F.rootInfo.bracket;
  Values(k,:)=[t(k),norm(S.Gd.er),A.T,F.Tstar,S.Tcmd,xx(36),(R.rhoB-A.hB)/R.rhoB, ...
   rad2deg(A.alpha),rad2deg(F.alpha),rad2deg(S.Air.alpha),rad2deg(S.Air.beta),S.Air.kappa, ...
   rad2deg(dcm_distance(S.Cc,S.Cr)),rad2deg(dcm_distance(S.Cr,S.Air.CIB)),rad2deg(dcm_distance(S.Cc,S.Air.CIB)), ...
   norm(S.Gd.vc)-R.Vc_min,hypot(S.Gd.vc(1),S.Gd.vc(2))-R.Vhc_min,max([0;-S.Gov.phi]),normPhi,corr,intervened,lam, ...
   S.Tcmd<P.Tmin-1e-8 || S.Tcmd>P.Tmax+1e-8,abs(preTdot)>P.TRateMax+1e-8,any(abs(preDelta)>P.deltaRateMax+1e-8),any(S.Alloc.active),norm(S.dAlloc),norm(S.dServo)];
 end
 assert(max(abs(er(:)-B.Y.er(:)))<1e-9);
 key=[r.method '_' r.case];RMS=sqrt(mean(Values(:,2).^2));complete=abs(t(end)-r.duration_s)<1e-10;
 if t(end)>t(1),RMStime=sqrt(trapz(t,Values(:,2).^2)/(t(end)-t(1)));else,RMStime=NaN;end
 late=t>=B.C.(['tSteady' upper(r.case(1)) r.case(2:end)]);
 if any(late),lateRMS=sqrt(mean(Values(late,2).^2));else,lateRMS=NaN;end
 holdDiff=NaN;
 if strcmp(r.case,'hold') && strcmp(r.method,'Transparent') && complete
  Z=load(fullfile(fileparts(H),'revision_wp3_20261107','results','Proposed_hold_export.mat'),'x','t');
  assert(isequal(size(x),size(Z.x)) && max(abs(t-Z.t))<1e-12);holdDiff=max(abs(x(:)-Z.x(:)));assert(holdDiff<=O.HoldStateMatchTol);
 end
 summary=struct('Case',r.case,'Method',r.method,'Source',r.path,'Samples',n,'Completed',complete,'HorizonLabel','partial', ...
  'ExpectedDuration_s',r.duration_s,'LastReturnedTime_s',t(end),'StoredRho',storedRho,'AuditRho',.10,'InitialStateDifference',initDiff, ...
  'PositionRMS_grid_m',RMS,'PositionRMS_time_m',RMStime,'PositionPeak_m',max(Values(:,2)),'PositionRMS_late_m',lateRMS, ...
  'ReferenceNonbackupStateViolationSamples',sum(any(Ref(:,1:25)<-1e-8,2)),'ReferenceThrustReserveViolationSamples',sum(any(Ref(:,24:25)<-1e-8,2)), ...
  'ReferenceEndpointReserveViolationSamples',sum(any(Ref(:,20:21)<-1e-8,2)),'ReferenceRootFailures',sum(~referenceRoot),'FeedbackRootFailures',sum(~feedbackRoot), ...
  'RawDerivativeNonbackupViolationSamples',sum(any(Phi(:,1:17)<-1e-8,2)),'NormalizedDerivativeViolationSamples',sum(Values(:,19)>O.NormalizedPhiTol), ...
  'BackupOutsideSamples',sum(Ref(:,26)<-O.BackupTol),'BackupDerivativeViolationSamples',sum(Phi(:,18)<-1e-8), ...
  'FeedbackRawThrustViolationSamples',sum(any(Feed(:,12:13)<-1e-8,2)),'ModeledActuatorViolationSamples',sum(any(Phys<-1e-8,2)), ...
  'StudyDomainViolationSamples',sum(any(Domain<-1e-8,2)),'ThrustTargetClippedSamples',sum(Values(:,23)), ...
  'ThrustRateLimitedSamples',sum(Values(:,24)),'SurfaceRateLimitedSamples',sum(Values(:,25)),'AllocatorBoundActiveSamples',sum(Values(:,26)), ...
  'ReferenceThrustMin_N',min(Values(:,3)),'ReferenceThrustMax_N',max(Values(:,3)),'RawThrustMin_N',min(Values(:,5)),'RawThrustMax_N',max(Values(:,5)), ...
  'ActuatorThrustMin_N',min(Values(:,6)),'ActuatorThrustMax_N',max(Values(:,6)),'MaxVBoverRho',max(Values(:,7)), ...
  'MinFeedbackRawEndpointLower_N',min(Feed(:,4)),'MinFeedbackRawEndpointUpper_N',min(Feed(:,5)), ...
  'MaxS4AttitudeSum_deg',max(Values(:,13)+Values(:,14)),'MaxDirectAttitude_deg',max(Values(:,15)), ...
  'GovernorIntervenedSamples',sum(Values(:,21)),'MaxAppliedU0Difference',max(max(abs(uApplied-u0))), ...
  'QPFallbackSamples',sum(qpFallback),'MaxQPKKT',max(QPKKT,[],1),'HoldStateDifference',holdDiff,'NumericalSampleTolerance',1e-8);
 if complete,summary.HorizonLabel='full';end
 allSummary{j}=summary;
 phinames=A.names;phiunits={'m/s3','m/s3','rad/s2','rad/s2','m/s3','m/s3','rad/s2','rad/s2','rad/s2','rad/s2','m/s3','m/s3','m/s3','N/s','N/s','N/s','N/s','weighted_set/s'};
 allRows=[allRows;rows(r,'reference',Ref,refnames,refunits,t);rows(r,'reference_derivative',Phi,phinames,phiunits,t); ...
  rows(r,'feedback',Feed,feednames,feedunits,t);rows(r,'modeled_actuator',Phys,physnames,physunits,t); ...
  rows(r,'study_domain',Domain,domnames,{'m/s','m/s','rad','rad','rad','rad'},t)]; %#ok<AGROW>
 save(fullfile(out,[key '_layered.mat']),'t','Ref','Feed','Phi','Phys','Domain','Values','valueNames','refnames','feednames','physnames','domnames','u0','uApplied','QPKKT','qpaths','qpFallback','summary','-v7');
 writetable(array2table(Values,'VariableNames',valueNames),fullfile(out,[key '_layered_timeseries.csv']));
 fprintf('WP4_AUDIT %s n=%d RMS=%.10f ref_reserve=%d KB_out=%d rawT_out=%d\n',key,n,RMS,summary.ReferenceThrustReserveViolationSamples,summary.BackupOutsideSamples,summary.FeedbackRawThrustViolationSamples);
end
Summary=[allSummary{:}];Margins=[allRows{:}];
fid=fopen(fullfile(out,'case_summary.json'),'w');fprintf(fid,'%s',jsonencode(Summary));fclose(fid);
fid=fopen(fullfile(out,'margin_summary.json'),'w');fprintf(fid,'%s',jsonencode(Margins));fclose(fid);
writetable(struct2table(Summary),fullfile(out,'case_summary.csv'));writetable(struct2table(Margins),fullfile(out,'margin_summary.csv'));
fprintf('WP4_ALL_NINE_LAYERED_AUDIT_COMPLETE\n');
end

function out=rows(r,category,A,names,units,t)
out=cell(size(A,2),1);
for j=1:size(A,2)
 [minimum,k]=min(A(:,j));negative=A(:,j)<-1e-8;f=find(negative,1);if isempty(f),first=NaN;else,first=t(f);end
 out{j}=struct('Case',r.case,'Method',r.method,'Category',category,'Margin',names{j},'Unit',units{j}, ...
  'Minimum',minimum,'MinimumTime_s',t(k),'NegativeSamples',sum(negative),'FirstNegativeTime_s',first,'Tolerance',1e-8);
end
end
