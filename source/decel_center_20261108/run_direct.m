function run_direct(scenario,step)
assert(ismember(scenario,{'C1','C2'}) && ismember(step,[.005 .01]));
here=fileparts(mfilename('fullpath')); root=fileparts(here); prev=fullfile(root,'decel_confirm_20261108');
src=fullfile(root,'revision_probe_20261106','source_snapshot');
addpath(here,prev,fullfile(root,'decel_d0_20261108'),fullfile(root,'decel_region_20261108'),...
    fullfile(root,'contribution_probe_20261106','full6dof'),fullfile(root,'contribution_probe_20261106','aoa'),fullfile(src,'core'),fullfile(src,'W2D'));
tag=sprintf('%s_%03d',scenario,round(step*1000));
folder=fullfile(here,'runs',tag);if ~exist(folder,'dir'),mkdir(folder);end
% Reuse the paired F configuration and independently generate the trim state.
B=load(fullfile(prev,'runs',sprintf('%s_F_%03d',scenario,round(step*1000)),'F.mat'),'C','O','Ts','Xs');
C=B.C;O=B.O;
[x0,I0]=cl6_initial_state('hold',C);assert(isequal(x0,B.Xs(1,:)') && isequal(x0(20:25),direct_reference(0,C.R)));
assert(C.MaxStep==step && C.R.rhoB==.1 && O.Horizon==35 && O.SpeedRise==10 && O.Course==0);
assert(O.Vlow==31+.5*strcmp(scenario,'C2') && all(O.WindFun(0)==0) && all(O.WindFun(20)==0));
for tt=unique([linspace(0,35,141),3-1e-6,3+1e-6,25-1e-6,25+1e-6])
    [xf,uc]=direct_reference(tt,C.R);[~,ub,h,zc]=moving_backup_affine(tt,xf,C.R);
    assert(isequal(zc,xf([1 2 4 5 6])) && isequal(uc,ub) && h==.01);
end
% Regression at identical plant/filter/servo states as the archived controller.
for tt=[0 3 8 14 20 25 35]
    k=find(abs(B.Ts-tt)<1e-9,1);xx=B.Xs(k,:)';xx(20:25)=direct_reference(tt,C.R);
    d=direct_snapshot(tt,xx,'turn',C,O);f=d0_snapshot(tt,xx,'turn',C,'F',O);
    for field={'Cmd','Cc','Cr','zeDot','whatDot','qrDot','omegaRDot','deltaCmd','deltaDot','Tcmd','TDot'}
        key=field{1};assert(isequaln(d.(key),f.(key)),key);
    end
end
clear B;
keep=[1:19,26:36]; y0=x0(keep);
global CP6_EXTENSION CP6_EVAL_TIME;CP6_EXTENSION=struct('reference',0,'feedback',0,'last_time_s',NaN);
grid=(0:C.sample_dt:O.Horizon)';
opts=odeset('RelTol',C.RelTol,'AbsTol',C.AbsTol,'MaxStep',C.MaxStep,'Events',@events);
fprintf('CENTER_START %s\n',tag);wall=tic;
[Ts,Y,te,ye,ie]=ode45(@rhs,grid,y0,opts);elapsed=toc(wall);
term=find(ie<=4,1);
if isempty(term)
    assert(abs(Ts(end)-O.Horizon)<1e-9,'Incomplete integration without an event');
    status='completed';
else
    status='root_boundary';
    if abs(Ts(end)-te(term))>1e-9
        assert(Ts(end)<te(term));Ts(end+1)=te(term);Y(end+1,:)=ye(term,:);
    end
end
Xs=zeros(numel(Ts),36);
names={'time','vground','vrequest','vref_ground','airSpeed','refAirSpeed','positionError','refLo','refHi','fbLo','fbHi',...
       'Tref','Traw','Tact','alphaActual_deg','eV','egamma_deg','echi_deg','theta_sum_deg','kappa','chartV','chartVh',...
       'surfaceAmp','surfaceRate','thrustAmp','thrustRate','surfaceCmdAtBound','thrustCmdClip','normalMargin','bankMargin'};
D=zeros(numel(Ts),numel(names));P=C.P7;tol=1e-8;
for k=1:numel(Ts)
    t=Ts(k);x=fullstate(t,Y(k,:)');Xs(k,:)=x'; S=direct_snapshot(t,x,'turn',C,O);xf=x(20:25);
    h=d0_raw_margins(t,x,C,O);F=rg_feedforward(xf,C.R);E=S.Cmd.e;vv=S.vfI+x(17:19);delta=x(33:35);T=x(36);
    D(k,:)=[t,hypot(S.Air.vgI(1),S.Air.vgI(2)),norm(S.Gd.rdDot(1:2)),norm(vv(1:2)),S.Air.V,xf(1),...
         norm(S.Gd.er),h',F.T,S.Tcmd,T,rad2deg(S.Air.alpha),E(1),rad2deg(E(2:3))',rad2deg(dist(S.Cc,S.Cr)+dist(S.Cr,S.Air.CIB)),...
         S.Air.kappa,norm(S.Gd.vc),norm(S.Gd.vc(1:2)),any(delta<=P.deltaMin+tol|delta>=P.deltaMax-tol),...
         any(abs(S.deltaDot)>=P.deltaRateMax-tol),T<=P.Tmin+tol||T>=P.Tmax-tol,abs(S.TDot)>=P.TRateMax-tol,...
         any(S.deltaCmd<=P.deltaMin+tol|S.deltaCmd>=P.deltaMax-tol),S.Tcmd<0||S.Tcmd>60,S.Cmd.aGamma-C.R.aGammaMin,...
         tan(C.R.muMax)*S.Cmd.aGamma-abs(S.Cmd.aChi)];
end
Summary=struct('case',scenario,'step',step,'status',status,'end_time',Ts(end),'events_s',te,'event_indices',ie,...
    'samples',numel(Ts),'rhs_extension',CP6_EXTENSION,'wall_seconds',elapsed);
save(fullfile(folder,'direct.mat'),'Ts','Xs','D','names','C','O','I0','Summary','-v7');
writetable(array2table(D,'VariableNames',names),fullfile(folder,'trace.csv'));
fid=fopen(fullfile(folder,'summary.json'),'w');fprintf(fid,'%s',jsonencode(Summary));fclose(fid);
fprintf('CENTER_DONE %s %s %.12f\n',tag,status,Ts(end));
    function x=fullstate(t,y)
        x=[y(1:19);direct_reference(t,C.R);y(20:30)];
    end
    function dy=rhs(t,y)
        CP6_EVAL_TIME=t;x=fullstate(t,y);S=direct_snapshot(t,x,'turn',C,O);
        U=struct('delta',x(33:35),'T',x(36),'wind_fun',O.WindFun);
        dy=[t2_gna_6dof_quat(t,x(1:13),U,C.P);S.zeDot;S.whatDot;S.qrDot;S.omegaRDot;S.deltaDot;S.TDot];
    end
    function [g,stop,dir]=events(t,y)
        x=fullstate(t,y);G=cp6_guidance(t,x(1:3),x(14:16),x(17:19),C,O);
        g=[d0_raw_margins(t,x,C,O);norm(G.vc)-C.R.Vc_min;norm(G.vc(1:2))-C.R.Vhc_min];
        stop=[ones(4,1);0;0];dir=-ones(6,1);
    end
end
function a=dist(X,Y)
a=acos(min(1,max(-1,(trace(X'*Y)-1)/2)));
end
