function A=rev411_initial_admissibility_screen(Samples,R,C0)
%REV411_INITIAL_ADMISSIBILITY_SCREEN Check all randomized t=0 governor states.
%
% This is an assumption-consistency screen, not a performance filter.
% It does not reject uncertainty realizations: every sample is retained.
% The initialization itself is constructed to satisfy the governor's declared
% backup-feasible premise.

if nargin<3 || isempty(C0), C0=cl6_params(); end

N=numel(Samples);
minSigma=inf;
minPhi=inf;
minVraw=inf; maxVraw=-inf;
minVf=inf; maxVf=-inf;
nChart=0;
nBackup=0;

for k=1:N
    C=rev41_apply_sample(C0,Samples(k),R);
    [x0,I0]=cl6_initial_state(R.mode,C);

    r=x0(1:3);
    ze=x0(14:16);
    what=x0(17:19);
    xf=x0(20:25);

    Gd=cl6_guidance_command(0,r,ze,what,R.mode,C);
    vc=Gd.vc;

    vraw=norm(vc);
    vh=hypot(vc(1),vc(2));
    chart=(vraw>C.R.Vc_min) && (vh>C.R.Vhc_min);
    if chart, nChart=nChart+1; end

    if ~chart
        continue;
    end

    [cc,~]=rg_raw_command_chart(vc,C.R);
    u0=rg_nominal_control(xf,cc,C.R);
    Aff=rg_constraint_affine(xf,C.R);
    uB=rg_backup_control(xf,C.R);
    sigma=Aff.c+Aff.G*uB;

    backup=min(sigma)>=-1e-9;
    if backup, nBackup=nBackup+1; end

    minSigma=min(minSigma,min(sigma));

    if backup
        Gov=rg_scalar_governor(xf,u0,C.R);
        minPhi=min(minPhi,min(Gov.phi));
    end

    minVraw=min(minVraw,vraw);
    maxVraw=max(maxVraw,vraw);
    minVf=min(minVf,xf(1));
    maxVf=max(maxVf,xf(1));

    % Initialization audit: filter starts at declared backup trim.
    if abs(xf(1)-C.R.Vb)>1e-10 || abs(xf(2)-C.R.gammaB)>1e-10 || ...
            norm(xf(4:6))>1e-10
        error('Rev411:InitState','Trial %d did not start at the backup filter state.',k);
    end

    %#ok<NASGU>
    I0=I0;
end

fprintf('\nREV.4.1.1 INITIAL-ADMISSIBILITY SCREEN\n');
fprintf('======================================\n');
fprintf('raw-command chart regular = %d / %d\n',nChart,N);
fprintf('backup-feasible x_f(0)     = %d / %d\n',nBackup,N);
fprintf('min initial backup sigma   = %.9e\n',minSigma);
fprintf('min initial governed phi   = %.9e\n',minPhi);
fprintf('raw-command speed range    = [%.3f %.3f] m/s\n',minVraw,maxVraw);
fprintf('governed initial V_f range = [%.3f %.3f] m/s\n',minVf,maxVf);

if nChart~=N
    error('Rev411:RawChart','At least one raw t=0 command is outside the chart.');
end
if nBackup~=N || minSigma<-1e-9
    error('Rev411:Backup','At least one randomized trial is not initially backup-feasible.');
end

A=struct('N',N,'nChart',nChart,'nBackup',nBackup, ...
    'minSigma',minSigma,'minPhi',minPhi, ...
    'rawVRange',[minVraw maxVraw],'VfRange',[minVf maxVf], ...
    'status','REV411_INITIAL_ADMISSIBILITY_ALL_SAMPLES_PASSED');
end
