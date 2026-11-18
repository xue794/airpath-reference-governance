function Q=a3_governor(x,u0,R,arm,O)
if strcmp(arm,'N10')
    Q=w2e_governor(x,u0,R,'MatchedQP',O,[]);
    Q.appliedC=Q.C.c; Q.appliedG=Q.C.G; Q.angleH=[NaN;NaN];
    return
end
A=rg_constraint_affine(x,R); c=A.c; G=A.G; angleH=[NaN;NaN];
switch arm
    case 'D16'
        keep=[1:13,16:18]; c=c(keep); G=G(keep,:);
    case 'Am'
        angle=cp_aoa_rows(x,R,'Am',A.alpha);
        c(14:15)=angle.c; G(14:15,:)=angle.G; angleH=angle.h;
    otherwise
        error('A3:Arm','Unexpected arm %s',arm);
end
uB=rg_backup_control(x,R); opt=O.QP; opt.FeasibleStart=uB;
[u,solver]=w2b_qp3(u0,-G,c,opt.Scale,opt);
if ~solver.Success,error('A3:QPSolver','%s',solver.Status);end
Q=struct('u',u,'uB',uB,'C',A,'Solver',solver,'appliedC',c,'appliedG',G,...
    'phi',c+G*u,'sigma',c+G*uB,'angleH',angleH);
end
