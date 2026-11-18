function a3_run_one(h,job,O)
id=job.id; arm=job.arm; step=job.step; level='coarse';if step<.01,level='fine';end
key=sprintf('id_%03d_%s_%s',id,arm,level); base=fullfile(h,'runs',key);
try
    prior=fullfile(fileparts(h),'a1','runs',sprintf('id_%03d',id));
    L=load([prior '.mat'],'S'); B=L.S; C=B.C; mode=O.Mode; x0=B.x(1,:).'; grid=B.t(:);
    seen=1; ts=zeros(size(grid)); xs=zeros(numel(grid),36); xs(1,:)=x0.';
    attempt=0; attemptedState=x0; lastlog=-20; started=tic; eid='';emsg='';
    reused=strcmp(arm,'N10') && step==.01;
    if reused
        t=grid; x=B.x; status='archived_completed'; attempt=NaN; attemptedState=[]; wall=0;
    else
        opt=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',step,'OutputFcn',@monitor);
        try
            [t,x]=ode45(@rhs,grid,x0,opt);
            if t(end)==grid(end),status='completed';else,status='solver_stopped';end
        catch failed
            t=ts(1:seen);x=xs(1:seen,:);status='terminated';eid=failed.identifier;emsg=getReport(failed,'extended','hyperlinks','off');
        end
        wall=toc(started);
    end
    save([base '_states.mat'],'id','arm','level','step','t','x','status','attempt','attemptedState','eid','emsg','wall','-v7');
    names={'time_s','pos_error_m','ref_T_N','ref_raw_lower_N','ref_raw_upper_N','VB_over_rho',...
       'fb_raw_lower_N','fb_raw_upper_N','T_raw_N','T_act_N','T_rate_Nps','T_target_N','surface_margin_rad',...
       'surface_rate_margin_radps','T_rate_margin_Nps','violation','kkt','backup_violation','Am_backup_lo_radps',...
       'Am_backup_hi_radps','Am_h_lo_rad','Am_h_hi_rad','u1','u2','u3','kappa','ref_alpha_rad','fb_alpha_rad'};
    data=nan(numel(t),numel(names)); P=C.P7;
    for k=1:numel(t)
        sn=a3_snapshot(t(k),x(k,:).',mode,C,arm,O); gov=sn.Gov; A=gov.C; F=sn.Cmd; Q=gov.Solver;
        lo=rg_N_map(C.R.alphaLo,sn.Air.V,F.X,C.R);hi=rg_N_map(C.R.alphaHi,sn.Air.V,F.X,C.R);
        rowviol=max([0;-gov.phi./(1+abs(gov.appliedC)+abs(gov.appliedG)*abs(gov.u))]);
        bviol=max([0;-gov.sigma./(1+abs(gov.appliedC)+abs(gov.appliedG)*abs(gov.uB))]);
        bk=[NaN;NaN]; if strcmp(arm,'Am'),bk=gov.sigma(14:15);end
        kkt=max([Q.PrimalNormalized,Q.StationarityNormalized,Q.ComplementarityNormalized,Q.DualViolationNormalized]);
        d=x(k,33:35).'; surf=min([d-P.deltaMin;P.deltaMax-d]);
        data(k,:)=[t(k),norm(sn.Gd.er),A.T,A.hNminus+C.R.DeltaN,A.hNplus+C.R.DeltaN,...
          (C.R.rhoB-A.hB)/C.R.rhoB,F.N-lo.value,hi.value-F.N,sn.Tcmd,x(k,36),sn.TDot,...
          min(P.Tmax,max(P.Tmin,sn.Tcmd)),surf,min(P.deltaRateMax-abs(sn.deltaDot)),P.TRateMax-abs(sn.TDot),...
          rowviol,kkt,bviol,bk.',gov.angleH.',gov.u.',sn.Air.kappa,A.alpha,F.alpha];
    end
    save([base '_samples.mat'],'names','data','-v7');
    if reused
        old=load([prior '_rows.mat'],'gapRef','gapFeedback');
        assert(isequal(data(:,23:25),B.Y.uGov) && isequal(data(:,9),B.Y.Tcmd) && ...
            isequal(data(:,4:5),old.gapRef) && isequal(data(:,7:8),old.gapFeedback));
    end
    record=struct('id',id,'arm',arm,'max_step',step,'status',status,'n',numel(t),'last_saved',t(end),...
        'last_rhs',attempt,'error_id',eid,'error_message',emsg,'wall_s',wall,'version',version,'archived_replay',reused);
    write_json([base '_status.json'],record);
catch ME
    write_json([base '_error.json'],struct('error',getReport(ME,'extended','hyperlinks','off')));
end
    function dx=rhs(tt,xx)
        attempt=tt;attemptedState=xx;dx=a3_rhs(tt,xx,mode,C,arm,O);
    end
    function stop=monitor(tt,xx,flag)
        stop=0;if ~isempty(flag),return;end
        for ii=1:numel(tt)
            if tt(ii)>ts(seen)+1e-12,seen=seen+1;ts(seen)=tt(ii);xs(seen,:)=xx(:,ii).';end
        end
        if ~isempty(tt)&&tt(end)-lastlog>=14.999
            lastlog=tt(end);write_json([base '_progress.json'],struct('last_saved',lastlog,'wall_s',toc(started)));
        end
    end
end

function write_json(p,s)
f=fopen([p '.tmp'],'w');assert(f>=0);fwrite(f,jsonencode(s),'char');fclose(f);movefile([p '.tmp'],p,'f');
end
