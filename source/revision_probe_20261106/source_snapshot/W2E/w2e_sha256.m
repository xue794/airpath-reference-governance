function h=w2e_sha256(file)
%W2E_SHA256 SHA-256 of a file using Java available in MATLAB R2018a.
fid=fopen(file,'rb');
if fid<0,error('W2E:Hash','Cannot open %s',file);end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while true
    b=fread(fid,1024*1024,'*uint8');
    if isempty(b),break;end
    md.update(typecast(b(:),'int8'));
end
d=typecast(md.digest(),'uint8');
h=lower(reshape(dec2hex(d,2).',1,[]));
end
