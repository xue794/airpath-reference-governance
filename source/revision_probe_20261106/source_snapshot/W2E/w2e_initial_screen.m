function A=w2e_initial_screen(Samples,Rmc,Cbase,O,outdir)
%W2E_INITIAL_SCREEN All 500 frozen samples, rhoB=0.10, no rejection sampling.
N=numel(Samples);Rows=cell(N,1);allPass=true;minSigma=inf;minScalarPhi=inf;minQPPhi=inf;maxVBover=0;
for k=1:N
 C=rev41_apply_sample(Cbase,Samples(k),Rmc);[x0,~]=cl6_initial_state(O.Mode,C);r=x0(1:3);ze=x0(14:16);what=x0(17:19);xf=x0(20:25);
 Gd=cl6_guidance_command(0,r,ze,what,O.Mode,C);vc=Gd.vc;vraw=norm(vc);vh=hypot(vc(1),vc(2));chart=(vraw>C.R.Vc_min)&&(vh>C.R.Vhc_min);
 VB=NaN;backupState=false;backupCandidate=false;scalarOK=false;qpOK=false;mnSig=NaN;mnP=NaN;mnQ=NaN;
 if chart
  [cc,~]=rg_raw_command_chart(vc,C.R);u0=rg_nominal_control(xf,cc,C.R);Aff=rg_constraint_affine(xf,C.R);uB=rg_backup_control(xf,C.R);sigma=Aff.c+Aff.G*uB;mnSig=min(sigma);xi=rg_backup_state(xf,C.R);VB=xi.'*C.R.PB*xi;backupState=VB<=C.R.rhoB+1e-12;backupCandidate=mnSig>=-1e-9;
  if backupState&&backupCandidate
   try,GP=rg_scalar_governor(xf,u0,C.R);scalarOK=true;mnP=min(GP.phi);catch,scalarOK=false;end
   try,GQ=rg_matched_qp(xf,u0,C.R,O.QP,Aff);qpOK=GQ.Success;mnQ=min(GQ.phi);catch,qpOK=false;end
  end
 end
 pass=chart&&backupState&&backupCandidate&&scalarOK&&qpOK;allPass=allPass&&pass;if isfinite(mnSig),minSigma=min(minSigma,mnSig);end;if isfinite(mnP),minScalarPhi=min(minScalarPhi,mnP);end;if isfinite(mnQ),minQPPhi=min(minQPPhi,mnQ);end;if isfinite(VB),maxVBover=max(maxVBover,VB/C.R.rhoB);end
 Rows{k}=struct('Index',k,'RawChartRegular',chart,'RawCommandSpeed_mps',vraw,'InitialVB_over_rhoB',VB/C.R.rhoB,'BackupStatePass',backupState,'BackupCandidatePass',backupCandidate,'MinBackupSigma',mnSig,'ProposedInitialSolvePass',scalarOK,'ProposedMinPhi',mnP,'QPInitialSolvePass',qpOK,'QPMinPhi',mnQ,'AllPass',pass);
end
writetable(struct2table([Rows{:}]),fullfile(outdir,'W2E_initial_admissibility.csv'));
A=struct('N',N,'AllPass',allPass,'MinBackupSigma',minSigma,'MinProposedPhi',minScalarPhi,'MinQPPhi',minQPPhi,'MaxInitialVB_over_rhoB',maxVBover,'Rows',{Rows});
if ~allPass,error('W2E:InitialScreen','At least one frozen sample failed the rhoB=0.10 initial admissibility/selector screen. No trial was rejected; simulation is halted.');end
end
