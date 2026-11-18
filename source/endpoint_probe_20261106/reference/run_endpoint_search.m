function run_endpoint_search()
here=fileparts(mfilename('fullpath'));task=fileparts(here);workspace=fileparts(task);
base=fullfile(workspace,'revision_probe_20261106','source_snapshot');addpath(fullfile(base,'core'));addpath(fullfile(base,'W2E'));addpath(here);
R=rg_params();R.rhoB=.10;O=w2b_options();rng(20261106);
for k=1:20
 z=randn(5,1);z=z/sqrt(z'*R.PB*z)*sqrt(.09)*rand();x=[R.Vb+z(1);z(2);0;z(3:5)];
 C=rg_constraint_affine(x,R);E=endpoint_affine_extended(x,R);
 assert(isequal(C.c,E.c)&&isequal(C.G,E.G),'Extended kernel differs inside original domain');
end
speeds=[30 31 32 34];courses=[30 60 90];results=struct([]);runs={};
for v=speeds
 for c=courses
  S=struct('speed',v,'gamma',0,'course',c);Q=sample_run(S,R,O,.025);
  results=[results;Q.summary];runs{end+1}=Q; %#ok<AGROW>
  fprintf('SCREEN V%g chi%g %s t=%.6f Tlast=%.6f eventT=%.6f\n',v,c,Q.summary.status,Q.t(end),Q.summary.lastT,Q.summary.eventThrust);
 end
end
save(fullfile(here,'screening.mat'),'results','runs','R','O');
write_json(fullfile(here,'screening_summary.json'),results);
writetable(struct2table(results),fullfile(here,'screening_summary.csv'));
fprintf('SCREENING_COMPLETE\n');
end
function Q=sample_run(S,R,O,dt)
x=[R.Vb;0;0;0;0;0];t=0;X=x';T=0;D=diagrow(x,R);U=[];status='complete';eventT=NaN;eventState=[];lastValidThrust=D(1,1);maxKKT=0;lastEventIndex=0;
for k=1:round(35/dt)
 d=endpoint_diagnostic(x,R);
 if ~d.root,status='root_lost_at_sample';break;end
 try,[u,sol]=endpoint_qp(t,x,S,R,O,'thrust15',false);catch ME,status=ME.identifier;break;end
 maxKKT=max(maxKKT,max([sol.PrimalNormalized,sol.StationarityNormalized,sol.ComplementarityNormalized]));
 xn=held(x,u,dt);dn=endpoint_diagnostic(xn,R);
 if min(dn.raw)<0
  % Exact held-input interpolation event, not a continuous-feedback event.
  lo=0;hi=dt;
  for it=1:45
   mid=(lo+hi)/2;dm=endpoint_diagnostic(held(x,u,mid),R);
   if min(dm.raw)>=0,lo=mid;else,hi=mid;end
  end
  tau=(lo+hi)/2;xn=held(x,u,tau);de=endpoint_diagnostic(xn,R);
  [~,id]=min(de.raw);eventT=de.endpointThrust(id);eventState=xn;lastEventIndex=id;
  t=t+tau;T=[T;t];X=[X;xn'];D=[D;diagrow(xn,R)];U=[U;u'];status='held_input_endpoint_event';break;
 end
 t=t+dt;x=xn;T=[T;t];X=[X;x'];D=[D;diagrow(x,R)];U=[U;u'];lastValidThrust=dn.T;
end
finiteT=D(isfinite(D(:,1)),1);
sumry=struct('speed',S.speed,'gamma',S.gamma,'course',S.course,'dt',dt,'status',status,'endTime',T(end), ...
'lastT',lastValidThrust,'eventThrust',eventT,'eventEndpointIndex',lastEventIndex,'minT',min(finiteT),'maxT',max(finiteT), ...
'thrustReserveViolations',sum(finiteT<3-1e-7|finiteT>57+1e-7),'conventionalViolations',sum(D(:,7)==0), ...
'minRawLower',min(D(:,2)),'minRawUpper',min(D(:,3)),'maxKKT',maxKKT);
Q=struct('t',T,'x',X,'u',U,'diagnostic',D,'eventState',eventState,'summary',sumry,'request',S);
end
function xn=held(x,u,dt),xn=[x(1:3)+dt*x(4:6)+.5*dt^2*u;x(4:6)+dt*u];end
function z=diagrow(x,R),d=endpoint_diagnostic(x,R);z=[d.T,d.raw',d.reserve',d.alpha,d.conventionalOK,min(d.psi)];end
function write_json(p,a),f=fopen(p,'w');fprintf(f,'%s',jsonencode(a));fclose(f);end
