function xi = rg_backup_state(x,R)
%RG_BACKUP_STATE Course-independent backup coordinates.
xi=[x(1)-R.Vb; x(2)-R.gammaB; x(4); x(5); x(6)];
end
