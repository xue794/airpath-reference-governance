function replay_confirm(scenario,method,step)
% Re-evaluate saved states. This routine has no time integration.
here=fileparts(fileparts(mfilename('fullpath')));root=fileparts(here);
addpath(here,fullfile(root,'decel_d0_20261108'),fullfile(root,'decel_region_20261108'),fullfile(root,'contribution_probe_20261106','aoa'),...
    fullfile(root,'contribution_probe_20261106','full6dof'),fullfile(root,'revision_probe_20261106','source_snapshot','core'),...
    fullfile(root,'revision_probe_20261106','source_snapshot','W2D'));
folder=fullfile(here,'runs',sprintf('%s_%s_%03d',scenario,method,round(step*1000)));
v=load(fullfile(folder,[method '.mat']));C=v.C;O=v.O;scale=O.QP.Scale(:);
names={'time','vground','vrequest','vref_ground','eV','egamma_deg','echi_deg','theta_sum_deg','kappa',...
    'hmove','u0V','u0G','u0C','uV','uG','uC','change','dual18','support18','support14','support15',...
    'primal','stationarity','dualViolation','complementarity','backupMinApplied',...
    'refLo','refHi','fbLo','fbHi','Tref','Traw','Tact','alphaActual_deg','positionError',...
    'surfaceAmp','surfaceRate','thrustAmp','thrustRate','surfaceCmdAtBound','thrustCmdClip',...
    'chartV','chartVh','normalMargin','bankMargin',...
    'D_scaledDelta','A_scaledDelta','deletedPhi14','deletedPhi15','FDKKT','FAKKT','airSpeed','refAirSpeed','windN','windE','windD','windDotN','windDotE','windDotD','whatN','whatE','whatD'};
X=nan(numel(v.Ts),numel(names));clear d0_governor;
for j=1:numel(v.Ts)
    t=v.Ts(j);xx=v.Xs(j,:)';xf=xx(20:25);
    S=d0_snapshot(t,xx,'turn',C,method,O);G=S.Gov;
    [~,~,hm]=moving_backup_affine(t,xf,C.R);
    mar=d0_raw_margins(t,xx,C,O);E=S.Cmd.e;sv=G.Solver;rows=G.rows(sv.SupportRows);
    lambda=0;if any(G.rows==18),lambda=sv.Dual(G.rows==18);end
    P=C.P7;delta=xx(33:35);T=xx(36);tol=1e-8;
    surfAmp=any(delta<=P.deltaMin+tol | delta>=P.deltaMax-tol);
    surfRate=any(abs(S.deltaDot)>=P.deltaRateMax-tol);
    dcActive=any(S.deltaCmd<=P.deltaMin+tol|S.deltaCmd>=P.deltaMax-tol);
    dr=[NaN NaN NaN NaN NaN NaN];
    if strcmp(method,'F')
        d=d0_governor(t,xf,S.u0,C.R,'D',O);a=d0_governor(t,xf,S.u0,C.R,'A',O);
        dd=norm((G.u-d.u)./scale,inf);da=norm((G.u-a.u)./scale,inf);
        dm=G.C.c(14:15)+G.C.G(14:15,:)*d.u;
        dr=[dd da dm' kkt(d.Solver) kkt(a.Solver)];
    end
    vv=S.vfI+xx(17:19);
    X(j,:)=[t,hypot(S.Air.vgI(1),S.Air.vgI(2)),hypot(S.Gd.rdDot(1),S.Gd.rdDot(2)),hypot(vv(1),vv(2)),...
        E(1),rad2deg(E(2:3))',rad2deg(dist(S.Cc,S.Cr)+dist(S.Cr,S.Air.CIB)),S.Air.kappa,hm,...
        S.u0',G.u',norm((G.u-S.u0)./scale,inf),lambda,any(rows==18),any(rows==14),any(rows==15),...
        sv.PrimalNormalized,sv.StationarityNormalized,sv.DualViolationNormalized,sv.ComplementarityNormalized,...
        min(G.cc+G.GG*G.backup),mar',G.C.T,S.Tcmd,T,rad2deg(S.Air.alpha),norm(S.Gd.er),...
        surfAmp,surfRate,(T<=P.Tmin+tol||T>=P.Tmax-tol),abs(S.TDot)>=P.TRateMax-tol,dcActive,(S.Tcmd<0||S.Tcmd>60),...
        S.chart.V,S.chart.Vh,S.Cmd.aGamma-C.R.aGammaMin,tan(C.R.muMax)*S.Cmd.aGamma-abs(S.Cmd.aChi),dr,S.Air.V,xf(1),S.Air.wI',S.Air.wDotI',xx(17:19)'];
end
writetable(array2table(X,'VariableNames',names),fullfile(folder,'replay.csv'));
fprintf('CONFIRM_REPLAY_OK %s %g %d\n',method,step,size(X,1));
end
function r=kkt(s)
r=max([s.PrimalNormalized s.StationarityNormalized s.DualViolationNormalized s.ComplementarityNormalized]);
end
function t=dist(a,b)
t=acos(min(1,max(-1,(trace(a'*b)-1)/2)));
end
