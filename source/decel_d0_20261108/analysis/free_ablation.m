function free_ablation()
% Same-state controls for the three arms with the backup row removed.
h=fileparts(fileparts(mfilename('fullpath')));r=fileparts(h);
addpath(h,fullfile(r,'contribution_probe_20261106','aoa'),fullfile(r,'contribution_probe_20261106','full6dof'),...
    fullfile(r,'revision_probe_20261106','source_snapshot','core'),fullfile(r,'revision_probe_20261106','source_snapshot','W2D'));
for id={'005','010'}
 dirp=fullfile(h,'runs',['N10_' id{1}]); v=load(fullfile(dirp,'N10.mat'));
 out=zeros(numel(v.Ts),7);opt=v.O;scale=opt.QP.Scale(:);
 for j=1:numel(v.Ts)
    t=v.Ts(j);x=v.Xs(j,:)';xf=x(20:25);
    S=d0_snapshot(t,x,'turn',v.C,'N10',opt);N=S.Gov;
    T=d0_governor(t,xf,S.u0,v.C.R,'T',opt);A=d0_governor(t,xf,S.u0,v.C.R,'Am',opt);
    row=N.C.c(14:15)+N.C.G(14:15,:)*T.u;
    out(j,:)=[t,norm((N.u-T.u)./scale,inf),norm((N.u-A.u)./scale,inf),row',kkt(T.Solver),kkt(A.Solver)];
 end
 writetable(array2table(out,'VariableNames',{'time','NT_delta','NA_delta','removed_phi14','removed_phi15','T_KKT','Am_KKT'}),fullfile(dirp,'free_replay.csv'));
end
fprintf('FREE_ABLATION_COMPLETE\n');
end
function v=kkt(s)
v=max([s.PrimalNormalized s.StationarityNormalized s.DualViolationNormalized s.ComplementarityNormalized]);
end
