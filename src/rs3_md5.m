function h = rs3_md5(str)
%RS3_MD5  Compute MD5 hex digest of a character vector.
% Uses Java MessageDigest (standard MATLAB).

if isstring(str)
    str = char(str);
end
if ~ischar(str)
    error('rs3_md5:InputMustBeChar', 'Input must be char or string.');
end

md = java.security.MessageDigest.getInstance('MD5');
md.update(uint8(str(:)'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d,2).',1,[]));
end
