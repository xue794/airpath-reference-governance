function audit_endpoint_validation()
here=fileparts(mfilename('fullpath'));task=fileparts(here);workspace=fileparts(task);
base=fullfile(workspace,'revision_probe_20261106','source_snapshot');addpath(fullfile(base,'core'));addpath(fullfile(base,'W2E'));addpath(here);
L=load(fullfile(here,'continuous_validation.mat'));E=struct([]);
for k=1:numel(L.Runs)
 Q=L.Runs{k};s=Q.summary;
 a=struct('method',Q.method,'maxStep',Q.maxStep,'eventRawMarginDerivative',NaN,'eventThrust',NaN, ...
 'eventSpeed',NaN,'eventGammaDeg',NaN,'eventCourseDeg',NaN,'eventCourseRateDegps',NaN,'eventBankDeg',NaN,'minimumMapDerivativeCoefficientLowerBound',Inf);
 for j=1:size(Q.x,1),a.minimumMapDerivativeCoefficientLowerBound=min(a.minimumMapDerivativeCoefficientLowerBound,derivative_lower(Q.x(j,:)',L.R));end
 if strcmp(s.status,'original_inversion_endpoint_event')
  x=Q.x(end,:)';[u,~,C]=endpoint_qp(Q.t(end),x,Q.request,L.R,L.O,Q.method,true);
  d=endpoint_diagnostic(x,L.R);row=13+s.eventEndpoint;deltaN=10;if strcmp(Q.method,'full17_N0'),deltaN=0;end
  a.eventRawMarginDerivative=C.c(row)+C.G(row,:)*u-L.R.kN*(d.raw(s.eventEndpoint)-deltaN);
  a.eventThrust=d.endpointThrust(s.eventEndpoint);a.eventSpeed=x(1);a.eventGammaDeg=x(2)*180/pi;a.eventCourseDeg=x(3)*180/pi;
  a.eventCourseRateDegps=x(6)*180/pi;a.eventBankDeg=d.bank*180/pi;
 end
 E=[E;a]; %#ok<AGROW>
end
V=struct([]);
for m={'thrust15','full17_N0','full17_N10'}
 ids=find(cellfun(@(q)strcmp(q.method,m{1}),L.Runs));a=L.Runs{ids(1)};b=L.Runs{ids(2)};
 [tc,ia,ib]=intersect(a.t,b.t); %#ok<ASGLU>
 difference=abs(a.x(ia,:)-b.x(ib,:));
 v=struct('method',m{1},'commonExactTimes',numel(tc),'maxStateDifferenceByComponent',max(difference,[],1), ...
 'endTimeDifference',abs(a.t(end)-b.t(end)),'minTDifference',abs(a.summary.minT-b.summary.minT), ...
 'maxTDifference',abs(a.summary.maxT-b.summary.maxT));V=[V;v]; %#ok<AGROW>
end
report=struct('eventAudits',E,'refinement',V,'scope','Numerical convergence and local outward derivative, not a formal safety proof.');
f=fopen(fullfile(here,'continuous_audit.json'),'w');fprintf(f,'%s',jsonencode(report));fclose(f);disp(report);fprintf('CONTINUOUS_AUDIT_COMPLETE\n');
end
function b=derivative_lower(x,R)
th=R.P.theta;hi=R.alphaHi;qS=.5*R.rho*R.P.S*x(1)^2;
assert(R.alphaLo==0 && hi>0 && hi<pi/2 && qS>0);
cld=pl([th(17),2*th(21),3*th(22),4*th(23)],hi);
cdd=pl([th(2),2*th(5),3*th(8),4*th(10)],hi);
cd=pl([th(1),th(2),th(5),th(8),th(10)],hi);
dx=qS*cd+R.P.m*(x(4)+R.g*sin(x(2)));
b=qS*cld+min(0,qS*cdd)*tan(hi)+min(dx,dx/cos(hi)^2);
end
function b=pl(c,h),b=c(1);for i=2:numel(c),b=b+min(0,c(i)*h^(i-1));end,end
