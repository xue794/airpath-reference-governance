function run_ep6_warm(Vlow,courseDegrees,tag,maxStep,turnRise,methods)
if nargin<1,Vlow=31;end
if nargin<2,courseDegrees=60;end
if nargin<3,tag=sprintf('V%g_C%g',Vlow,courseDegrees);end
if nargin<4,maxStep=.01;end
if nargin<5,turnRise=3;end
if nargin<6,methods={'T','N0','N10'};end
clear ep6_governor_warm;
here=fileparts(mfilename('fullpath'));workspace=fileparts(fileparts(here));
source=fullfile(workspace,'revision_probe_20261106','source_snapshot');
addpath(fullfile(source,'core'));addpath(fullfile(source,'W2D'));addpath(here);
outdir=fullfile(here,tag);if ~exist(outdir,'dir'),mkdir(outdir);end
C=cl6_params();assert(C.R.rhoB==.2);C.R.rhoB=.1;C.MaxStep=maxStep;
O=struct('QP',w2b_options(),'Vlow',Vlow,'SpeedStart',3,'SpeedRise',8,'TurnStart',16,'TurnRise',turnRise,'Course',courseDegrees*pi/180,'Horizon',35);
[O.Nodes,O.Weights]=gauss(24);[n48,w48]=gauss(48);O48=O;O48.Nodes=n48;O48.Weights=w48;
maxQuadError=0;maxDerivativeError=0;
for tt=linspace(0,O.Horizon,101)
 [p,v]=ep6_reference(tt,C,O);p48=ep6_reference(tt,C,O48);maxQuadError=max(maxQuadError,norm(p-p48,inf));
 if tt>1e-4,pPlus=ep6_reference(tt+1e-5,C,O);pMinus=ep6_reference(tt-1e-5,C,O);maxDerivativeError=max(maxDerivativeError,norm((pPlus-pMinus)/2e-5-v,inf));end
end
assert(maxQuadError<1e-9&&maxDerivativeError<1e-6,'EP6:Reference','Ground reference consistency check failed');
[x0,I0]=cl6_initial_state('hold',C);assert(all(ep6_raw_margins(0,x0,C)>0));
% Exact regular-domain kernel identity; extension only changes trial points
% beyond original inversion endpoints, where terminal event is being located.
for v=[31,35,38.1,42]
 xf=[v;0;0;0;0;0];A=rg_constraint_affine(xf,C.R);B=ep6_affine(xf,C.R);
 assert(isequal(A.c,B.c)&&isequal(A.G,B.G));
end
Results=cell(numel(methods),1);
for j=1:numel(methods)
 method=methods{j};Ts=[];Xs=[];lastTime=NaN;lastState=[];
 global EP6_EXTENSION EP6_EVAL_TIME;
 EP6_EXTENSION=struct('reference',0,'feedback',0,'last_time_s',NaN);EP6_EVAL_TIME=NaN;
 status='completed';errorId='';errorMessage='';eventTimes=[];eventStates=[];eventIndices=[];
 fprintf('EP6_START %s %s\n',tag,method);wall=tic;
 try
   grid=(0:C.sample_dt:O.Horizon)';
   opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep,'Events',@events,'OutputFcn',@capture);
   [t,x,te,xe,ie]=ode45(@rhs,grid,x0,opts);Ts=t;Xs=x;eventTimes=te;eventStates=xe;eventIndices=ie;
   if ~isempty(te),status='terminal_original_root_boundary';end
 catch ME
   status='interface_error_at_solver_trial';errorId=ME.identifier;errorMessage=ME.message;
 end
 wallSeconds=toc(wall);Extension=EP6_EXTENSION;
 [Ts,index]=unique(Ts,'stable');Xs=Xs(index,:);
 if ~isempty(eventTimes)&&abs(Ts(end)-eventTimes(end))>1e-9,Ts=[Ts;eventTimes(end)];Xs=[Xs;eventStates(end,:)];end
 N=numel(Ts);D=nan(N,29);
 for k=1:N
   tt=Ts(k);xx=Xs(k,:)';xf=xx(20:25);S=ep6_snapshot_warm(tt,xx,'turn',C,method,O);
   h=ep6_raw_margins(tt,xx,C);F=rg_feedforward(xf,C.R);
   if ~F.rootInfo.bracket
     % Boundary output may be outside by roundoff. Only evaluate its limiting
     % Tff using the same event continuation; flag original bracket separately.
     X=C.R.P.m*(xf(4)+C.R.g*sin(xf(2)));Nf=C.R.P.m*hypot(xf(1)*xf(5)+C.R.g*cos(xf(2)),xf(1)*cos(xf(2))*xf(6));
     al=ep6_alpha_root(xf(1),X,Nf,C.R,'reference');Aa=rg_nominal_aero(xf(1),al,C.R);Tf=(Aa.D+X)/cos(al);
   else,al=F.alpha;Tf=F.T;end
   xi=rg_backup_state(xf,C.R);psi=[xf(4)+C.R.k1*(xf(1)-C.R.Vmin);-xf(4)+C.R.k1*(C.R.Vmax-xf(1));xf(5)+C.R.k1*(xf(2)+C.R.gammaMax);-xf(5)+C.R.k1*(C.R.gammaMax-xf(2))];
   D(k,:)=[tt,norm(S.Gd.er),xf(1),S.Air.V,rad2deg(xf(3)),rad2deg(S.Air.chi),rad2deg(al),rad2deg(S.Cmd.alpha),rad2deg(S.Air.alpha), ...
       Tf,S.Tcmd,xx(36),h',xi'*C.R.PB*xi/C.R.rhoB,min(psi), ...
       C.R.Vmin<=xf(1)+1e-8&&xf(1)<=C.R.Vmax+1e-8,abs(xf(2))<=C.R.gammaMax+1e-8,all(abs(xf(4:6))<=C.R.nuMax+1e-8), ...
       F.aGamma>=C.R.aGammaMin-1e-8,abs(F.mu)<=C.R.muMax+1e-8,S.Gov.Solver.PrimalNormalized, ...
       F.rootInfo.bracket,S.Cmd.rootInfo.originalBracket,rad2deg(F.mu),rad2deg(xf(6)),rad2deg(S.Air.beta)];
 end
 names={'time_s','position_error_m','Vref_mps','Vactual_mps','course_ref_deg','course_actual_deg','alpha_ref_deg','alpha_command_deg','alpha_actual_deg', ...
     'Tref_N','Traw_command_N','Tactual_N','ref_raw_lower_N','ref_raw_upper_N','feedback_raw_lower_N','feedback_raw_upper_N','VB_over_rho','min_psi1', ...
     'speed_ok','climb_ok','rate_ok','normal_ok','bank_ok','QP_primal_normalized','original_reference_bracket','original_feedback_bracket','bank_ref_deg','course_rate_ref_degps','beta_actual_deg'};
 writetable(array2table(D,'VariableNames',names),fullfile(outdir,[method '_samples.csv']));
 failedTrialTime=NaN;if strcmp(status,'interface_error_at_solver_trial'),failedTrialTime=lastTime;end
 Summary=struct('method',method,'numerical_start','previous support proposal; unchanged w2b_qp3 verifies KKT','status',status,'error_id',errorId,'error_message',errorMessage,'end_time_s',Ts(end),'requested_time_s',O.Horizon, ...
   'event_times_s',eventTimes,'event_indices',eventIndices,'failed_trial_time_s',failedTrialTime,'last_rhs_evaluation_time_s',lastTime,'RMS_until_stop_m',sqrt(mean(D(:,2).^2)), ...
   'min_Tref_N',min(D(:,10)),'max_Tref_N',max(D(:,10)),'min_Traw_N',min(D(:,11)),'max_Traw_N',max(D(:,11)), ...
   'min_ref_raw_upper_N',min(D(:,14)),'min_feedback_raw_upper_N',min(D(:,16)),'max_alpha_actual_deg',max(D(:,9)), ...
   'conventional_violation_samples',sum(any(D(:,19:23)==0,2)),'min_psi1',min(D(:,18)), ...
   'reference_thrust_reserve_violation_samples',sum(D(:,10)<3-1e-6|D(:,10)>57+1e-6), ...
   'raw_thrust_hardware_violation_samples',sum(D(:,11)<-1e-6|D(:,11)>60+1e-6), ...
   'actual_thrust_hardware_violation_samples',sum(D(:,12)<-1e-6|D(:,12)>60+1e-6), ...
   'rootfree_margin_violation_samples',sum(any(D(:,13:16)<-1e-6,2)), ...
   'max_VB_over_rho',max(D(:,17)),'max_QP_primal_normalized',max(D(:,24)),'extension_evaluations',Extension, ...
   'wall_seconds',wallSeconds,'reference_quadrature_check_m',maxQuadError,'reference_derivative_check_mps',maxDerivativeError);
 Results{j}=Summary;save(fullfile(outdir,[method '.mat']),'Ts','Xs','D','Summary','C','O','I0','lastState','eventStates','-v7');
 fid=fopen(fullfile(outdir,[method '_summary.json']),'w');fprintf(fid,'%s',jsonencode(Summary));fclose(fid);disp(Summary);
end
fid=fopen(fullfile(outdir,'all_summary.json'),'w');fprintf(fid,'%s',jsonencode([Results{:}]));fclose(fid);fprintf('EP6_DONE\n');
    function dx=rhs(tt,xx)
       lastTime=tt;lastState=xx;EP6_EVAL_TIME=tt;
       S=ep6_snapshot_warm(tt,xx,'turn',C,method,O);U=cl6_plant_input(xx(33:35),xx(36),'hold',C);
       dp=t2_gna_6dof_quat(tt,xx(1:13),U,C.P);dx=[dp;S.zeDot;S.whatDot;S.xfDot;S.qrDot;S.omegaRDot;S.deltaDot;S.TDot];
    end
    function [v,term,dir]=events(tt,xx)
       v=ep6_raw_margins(tt,xx,C);term=ones(4,1);dir=-ones(4,1);
    end
    function stop=capture(tt,xx,flag)
       stop=0;if isempty(flag),Ts=[Ts;tt(:)];Xs=[Xs;xx'];end
    end
end
function [x,w]=gauss(n)
i=(1:n-1)';a=i./sqrt(4*i.^2-1);[V,D]=eig(diag(a,1)+diag(a,-1));[x,ix]=sort(diag(D));w=2*V(1,ix)'.^2;
end
