function chart_events()
% Resolve the first chart threshold from the last retained common-grid state.
% Extending the chart only for solver trial points localizes the original event.
here=fileparts(mfilename('fullpath'));study=fileparts(here);root=fileparts(study);
addpath(here,study,fullfile(root,'decel_region_20261108'),fullfile(root,'contribution_probe_20261106','aoa'),...
 fullfile(root,'contribution_probe_20261106','full6dof'),fullfile(root,'revision_probe_20261106','source_snapshot','core'),...
 fullfile(root,'revision_probe_20261106','source_snapshot','W2D'));
src=fileread(fullfile(study,'d0_snapshot.m'));
marker='[cc,chart]=rg_raw_command_chart(Gd.vc,C.R);';
assert(numel(strfind(src,marker))==1);
src=strrep(src,marker,'Rloc=C.R; Rloc.Vc_min=0.1; Rloc.Vhc_min=0.1; [cc,chart]=rg_raw_command_chart(Gd.vc,Rloc);');
src=strrep(src,'function S=d0_snapshot','function S=d0_chart_probe');
fid=fopen(fullfile(here,'d0_chart_probe.m'),'w');fprintf(fid,'%s',src);fclose(fid);rehash;
entries=dir(fullfile(study,'runs','*','*.mat'));
for i=1:numel(entries)
 L=load(fullfile(entries(i).folder,entries(i).name));
 if ~strcmp(L.Summary.error_id,'rg_raw_command_chart:singular'),continue;end
 C=L.C;O=L.O;arm=L.Summary.method;t0=L.Ts(end);x0=L.Xs(end,:)';
 runs=cell(2,1);
 for j=1:2
     step=[.001 .0005];opt=odeset('Events',@events,'MaxStep',step(j),'RelTol',1e-10,'AbsTol',1e-12);
     [~,~,te,xe,ie]=ode45(@rhs,[t0 t0+.12],x0,opt);
     assert(~isempty(te) && all(ie<=2));t=te(1);x=xe(1,:)';s=d0_chart_probe(t,x,'turn',C,arm,O);
     runs{j}=struct('time_s',t,'index',ie(1),'Vchart',s.chart.V,'Vhchart',s.chart.Vh,...
         'gaps',d0_raw_margins(t,x,C,O),'Tref',s.Gov.C.T,'Traw',s.Tcmd,'Tact',x(36));
 end
 result=struct('arm',arm,'last_saved',t0,'last_trial',L.Summary.last_rhs_evaluation_time_s,'runs',{runs});
 fid=fopen(fullfile(entries(i).folder,'chart_event.json'),'w');fprintf(fid,'%s',jsonencode(result));fclose(fid);
end
fprintf('CHART_EVENT_CHECK_DONE\n');
function dx=rhs(t,x)
 S=d0_chart_probe(t,x,'turn',C,arm,O);U=struct('delta',x(33:35),'T',x(36),'wind_fun',O.WindFun);
 dx=[t2_gna_6dof_quat(t,x(1:13),U,C.P);S.zeDot;S.whatDot;S.xfDot;S.qrDot;S.omegaRDot;S.deltaDot;S.TDot];
end
function [v,term,direction]=events(t,x)
 G=cp6_guidance(t,x(1:3),x(14:16),x(17:19),C,O);
 v=[norm(G.vc)-C.R.Vc_min;hypot(G.vc(1),G.vc(2))-C.R.Vhc_min;d0_raw_margins(t,x,C,O)];
 term=ones(6,1);direction=-term;
end
end
