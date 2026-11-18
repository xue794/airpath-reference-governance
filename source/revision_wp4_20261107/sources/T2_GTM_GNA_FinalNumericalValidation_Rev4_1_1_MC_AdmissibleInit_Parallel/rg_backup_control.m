function uB = rg_backup_control(x,R)
%RG_BACKUP_CONTROL Strictly feasible local backup feedback.
uB=[-R.kVp*(x(1)-R.Vb)-R.kVd*x(4);
    -R.kGammaP*(x(2)-R.gammaB)-R.kGammaD*x(5);
    -R.kChiD*x(6)];
end
