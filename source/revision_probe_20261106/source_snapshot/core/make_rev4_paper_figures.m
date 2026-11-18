function Files=make_rev4_paper_figures(Results4)
%MAKE_REV4_PAPER_FIGURES Publication-ready Rev.4 figures for MATLAB R2018a.
%
% Usage:
%   FilesPaper = make_rev4_paper_figures(Results4);
%
% The function does not rerun simulations or modify any numerical result.
% It only reformats the already accepted Rev.4.0.2 results for manuscript use.

if nargin<1 || isempty(Results4)
    error('Pass the accepted Results4 structure.');
end
if ~isfield(Results4,'status') || ...
        ~strcmp(Results4.status,'REV4_FINAL_NUMERICAL_VALIDATION_PASSED')
    error('Results4 is not an accepted Rev.4 result.');
end

root=fileparts(mfilename('fullpath'));
outdir=fullfile(root,'rev4_paper_figures');
if ~exist(outdir,'dir'), mkdir(outdir); end

V=Results4.Regression.CL.varying;
N=Results4.Regression.CL.noisy;
C=Results4.Certificate;

fs=9;
lw=1.15;
lw2=1.0;

Files={};

% -------------------------------------------------------------------------
% Fig. 1: position tracking, varying vs noisy.
% -------------------------------------------------------------------------
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.6]);

labs={'e_N [m]','e_E [m]','e_D [m]'};
for i=1:3
    subplot(3,1,i);
    plot(V.t,V.Y.er(:,i),'-','LineWidth',lw); hold on;
    plot(N.t,N.Y.er(:,i),'--','LineWidth',lw);
    grid on; box on;
    ylabel(labs{i},'Interpreter','none','FontSize',fs);
    set(gca,'FontSize',fs);
    xlim([min(V.t) max(V.t)]);
    if i==1
        legend({'varying wind','varying wind + measurement error'}, ...
            'Location','best','FontSize',fs-1);
    end
    if i==3
        xlabel('t [s]','Interpreter','none','FontSize',fs);
    end
end
Files{end+1}=local_export(f,outdir,'Fig_Rev4_PositionTracking'); %#ok<AGROW>
close(f);

% -------------------------------------------------------------------------
% Fig. 2: estimator + governor robustness.
% -------------------------------------------------------------------------
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.6]);

subplot(3,1,1);
plot(V.t,sqrt(sum(V.Y.tildew.^2,2)),'-','LineWidth',lw); hold on;
plot(N.t,sqrt(sum(N.Y.tildew.^2,2)),'--','LineWidth',lw);
grid on; box on;
ylabel('||w-w_hat|| [m/s]','Interpreter','none','FontSize',fs);
legend({'varying wind','with measurement error'},'Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

subplot(3,1,2);
plot(V.t,V.Y.lambda,'-','LineWidth',lw); hold on;
plot(N.t,N.Y.lambda,'--','LineWidth',lw);
grid on; box on;
ylabel('lambda*','Interpreter','none','FontSize',fs);
ylim([0 1.02]); set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

subplot(3,1,3);
plot(V.t,max(V.Y.minPhi,0),'-','LineWidth',lw); hold on;
plot(N.t,max(N.Y.minPhi,0),'--','LineWidth',lw);
grid on; box on;
ylabel('max(min phi_j,0)','Interpreter','none','FontSize',fs);
xlabel('t [s]','Interpreter','none','FontSize',fs);
set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

Files{end+1}=local_export(f,outdir,'Fig_Rev4_EstimatorGovernor');
close(f);

% -------------------------------------------------------------------------
% Fig. 3: physical inner-loop realization for varying wind.
% -------------------------------------------------------------------------
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 6.4]);

subplot(3,1,1);
plot(V.t,rad2deg(V.Y.thetaF),'-','LineWidth',lw); hold on;
plot(V.t,rad2deg(V.Y.thetaR),'--','LineWidth',lw);
plot(V.t,rad2deg(V.Y.thetaF+V.Y.thetaR),':','LineWidth',lw2);
grid on; box on;
ylabel('attitude error [deg]','Interpreter','none','FontSize',fs);
legend({'theta_f','theta_R','theta_f+theta_R'},'Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

subplot(3,1,2);
plot(V.t,rad2deg(V.Y.deltaAct),'LineWidth',lw);
grid on; box on;
ylabel('delta_act [deg]','Interpreter','none','FontSize',fs);
legend({'aileron','elevator','rudder'},'Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

subplot(3,1,3);
plot(V.t,V.Y.Tact,'-','LineWidth',lw); hold on;
plot(V.t,V.Y.Tcmd,'--','LineWidth',lw);
grid on; box on;
ylabel('thrust [N]','Interpreter','none','FontSize',fs);
xlabel('t [s]','Interpreter','none','FontSize',fs);
legend({'T_act','T_cmd'},'Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([min(V.t) max(V.t)]);

Files{end+1}=local_export(f,outdir,'Fig_Rev4_PhysicalRealization');
close(f);

% -------------------------------------------------------------------------
% Fig. 4: local macro certificate.
% Samples: 1=hold, 2:26=reference lifting, 27:45=realized trajectory.
% -------------------------------------------------------------------------
Nfam=numel(C.Family);
idx=(1:Nfam).';
maxR=zeros(Nfam,1);
q=zeros(Nfam,1);
for k=1:Nfam
    maxR(k)=C.Family(k).maxReal;
    q(k)=C.Family(k).q;
end

f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 5.2]);

subplot(2,1,1);
plot(idx,maxR,'o-','LineWidth',lw,'MarkerSize',4); hold on;
yl=ylim;
plot([26.5 26.5],yl,'--','LineWidth',0.8);
grid on; box on;
ylabel('max Re eig(A_C)','Interpreter','none','FontSize',fs);
set(gca,'FontSize',fs); xlim([1 Nfam]);
text(8,yl(2)-0.10*(yl(2)-yl(1)),'reference lifting','FontSize',fs-1);
text(31,yl(2)-0.10*(yl(2)-yl(1)),'realized 6DOF','FontSize',fs-1);

subplot(2,1,2);
plot(idx,q,'o-','LineWidth',lw,'MarkerSize',4); hold on;
yl=ylim;
plot([26.5 26.5],yl,'--','LineWidth',0.8);
grid on; box on;
ylabel('min eig(Q_C)','Interpreter','none','FontSize',fs);
xlabel('sample index','Interpreter','none','FontSize',fs);
set(gca,'FontSize',fs); xlim([1 Nfam]);

Files{end+1}=local_export(f,outdir,'Fig_Rev4_MacroCertificate');
close(f);

fprintf('\nRev.4 manuscript figures written to:\n  %s\n',outdir);
end

function F=local_export(fig,outdir,stem)
png=fullfile(outdir,[stem '.png']);
epsf=fullfile(outdir,[stem '.eps']);
print(fig,png,'-dpng','-r300');
print(fig,epsf,'-depsc2','-r600');
F=struct('png',png,'eps',epsf);
end
