function w2f_write_csv(filename,T)
%W2F_WRITE_CSV R2018a CSV output helper.
if istable(T)
    writetable(T,filename);
else
    error('W2F:CSV','Expected a table.');
end
end
