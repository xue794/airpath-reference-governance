function run_endpoint_continuous()
% Continuous-feedback validation after the explicitly exploratory finite screen.
% Same C2 request and same weighted QP. Only endpoint rows/reserve change.
here=fileparts(mfilename('fullpath'));task=fileparts(here);workspace=fileparts(task);
base=fullfile(workspace,'revision_probe_20261106','source_snapshot');addpath(fullfile(base,'core'));addpath(fullfile(base,'W2E'));addpath(here);
R=rg_params();R.rhoB=.10;O=w2b_options();S=struct('speed',31,'gamma',0,'course',60);
methods={'thrust15','full17_N0','full17_N10'};steps=[.01 .005];results=struct([]);Runs={};
for ds=steps
 for k=1:numel(methods)
  Q=one_run(S,methods{k},ds,R,O);results=[results;Q.summary];Runs{end+1}=Q; %#ok<AGROW>
  prefix=sprintf('continuous_%s_maxstep%g',methods{k},ds);
  names={'t','V','gamma','chi','nuV','nuGamma','nuChi','Tff','rawNlo','rawNhi','reserveNlo','reserveNhi','alpha','conventionalOK','minPsi1','uV','uGamma','uChi','endpointSupport','thrustSupport'};
  writetable(array2table([Q.t Q.x Q.diag Q.u Q.support],'VariableNames',names),fullfile(here,[prefix '.csv']));
  fprintf('CONTINUOUS %s maxstep=%g status=%s tend=%.9f Tmin=%.9f Tmax=%.9f minNhi=%.9g extension=%d\n',methods{k},ds,Q.summary.status,Q.t(end),Q.summary.minT,Q.summary.maxT,Q.summary.minRawNhi,Q.summary.extensionCalls);
 end
end
save(fullfile(here,'continuous_validation.mat'),'Runs','results','R','O','S');
f=fopen(fullfile(here,'continuous_summary.json'),'w');fprintf(f,'%s',jsonencode(results));fclose(f);
writetable(struct2table(results),fullfile(here,'continuous_summary.csv'));fprintf('CONTINUOUS_VALIDATION_COMPLETE\n');
end
function Q=one_run(S,method,ds,R,O)
extensionCalls=0;evalCalls=0;lastEvalTime=NaN;acceptedTimes=[];acceptedStates=[];
initial=[R.Vb;0;0;0;0;0];status='complete';err='';te=[];xe=[];ie=[];
opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',ds,'Events',@events,'OutputFcn',@capture);
try
 [t,x,te,xe,ie]=ode45(@rhs,(0:.02:35)',initial,opts);
 if ~isempty(te),status='original_inversion_endpoint_event';end
catch ME
 status='solver_interface_error';err=[ME.identifier ' ' ME.message];t=acceptedTimes;x=acceptedStates;
end
if ~isempty(te)&&abs(t(end)-te(end))>1e-10,t=[t;te(end)];x=[x;xe(end,:)];end
[t,ix]=unique(t,'stable');x=x(ix,:);n=numel(t);z=nan(n,8);u=nan(n,3);supports=zeros(n,2);maxKKT=0;rootFail=0;
for k=1:n
 d=endpoint_diagnostic(x(k,:)',R);Tff=d.T;
 if ~d.root
  rootFail=rootFail+1;
  [mm,ii]=min(abs(d.raw));if mm<1e-6,Tff=d.endpointThrust(ii);end
 end
 z(k,:)=[Tff,d.raw',d.reserve',d.alpha,d.conventionalOK,min(d.psi)];
 [uk,D]=endpoint_qp(t(k),x(k,:)',S,R,O,method,true);u(k,:)=uk';
 supports(k,:)=[any(ismember(D.supportOriginal,[14 15])),any(ismember(D.supportOriginal,[16 17]))];
 maxKKT=max(maxKKT,max([D.PrimalNormalized,D.StationarityNormalized,D.ComplementarityNormalized]));
end
if isempty(te),eventTime=NaN;eventThrust=NaN;eventIndex=0;eventOutwardRate=NaN;
else
 eventTime=te(end);eventIndex=ie(end);ed=endpoint_diagnostic(xe(end,:)',R);eventThrust=ed.endpointThrust(eventIndex);
 [ue,~,Ce]=endpoint_qp(eventTime,xe(end,:)',S,R,O,method,true);
 eventRow=13+eventIndex;eventDeltaN=R.DeltaN;if strcmp(method,'full17_N0'),eventDeltaN=0;end
 eventOutwardRate=Ce.c(eventRow)+Ce.G(eventRow,:)*ue-R.kN*(ed.raw(eventIndex)-eventDeltaN);
end
vErr=x(:,1)-S.speed;chiErr=x(:,3)-S.course*pi/180;
summary=struct('method',method,'maxStep',ds,'status',status,'error',err,'endTime',t(end),'samples',n, ...
 'eventTime',eventTime,'eventEndpoint',eventIndex,'eventThrust',eventThrust,'eventRawMarginDerivative',eventOutwardRate,'extensionCalls',extensionCalls,'RHSCalls',evalCalls,'lastTrialTime',lastEvalTime, ...
 'minT',min(z(:,1)),'maxT',max(z(:,1)),'thrustReserveViolations',sum(z(:,1)<3-1e-7|z(:,1)>57+1e-7), ...
 'minRawNlo',min(z(:,2)),'minRawNhi',min(z(:,3)),'minReserveNlo',min(z(:,4)),'minReserveNhi',min(z(:,5)), ...
 'originalBracketFalseSamples',rootFail,'originalDomainBeyondToleranceSamples',sum(any(z(:,2:3)<-1e-6,2)), ...
 'conventionalViolations',sum(z(:,7)==0),'minPsi1',min(z(:,8)),'maxKKT',maxKKT, ...
 'endpointSupportSamples',sum(supports(:,1)),'thrustSupportSamples',sum(supports(:,2)), ...
 'finalV',x(end,1),'finalGammaDeg',x(end,2)*180/pi,'finalCourseDeg',x(end,3)*180/pi, ...
 'post16SpeedTargetRMS',sqrt(mean(vErr(t>=16).^2)),'post16CourseTargetRMSdeg',sqrt(mean(chiErr(t>=16).^2))*180/pi);
Q=struct('request',S,'method',method,'maxStep',ds,'t',t,'x',x,'diag',z,'u',u,'support',supports,'summary',summary);
    function dx=rhs(tt,xx)
        evalCalls=evalCalls+1;lastEvalTime=tt;
        vv=raw_endpoints(xx,R);if min(vv)<0,extensionCalls=extensionCalls+1;end
        uu=endpoint_qp(tt,xx,S,R,O,method,true);dx=[xx(4:6);uu];
    end
    function [v,terminal,dir]=events(~,xx)
        v=raw_endpoints(xx,R);terminal=[1;1];dir=[-1;-1];
    end
    function stop=capture(tt,xx,flag)
        stop=0;if isempty(flag),acceptedTimes=[acceptedTimes;tt(:)];acceptedStates=[acceptedStates;xx'];end
    end
end
function v=raw_endpoints(x,R)
X=R.P.m*(x(4)+R.g*sin(x(2)));ag=x(1)*x(5)+R.g*cos(x(2));ac=x(1)*cos(x(2))*x(6);N=R.P.m*hypot(ag,ac);
lo=rg_N_map(R.alphaLo,x(1),X,R);hi=rg_N_map(R.alphaHi,x(1),X,R);v=[N-lo.value;hi.value-N];
end
