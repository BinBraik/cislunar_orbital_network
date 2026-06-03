function s = sanitize_fname(s)
%SANITIZE_FNAME  Make a string safe for filenames.
s = strrep(s,' ','_');
s = regexprep(s,'[^a-zA-Z0-9_\-]','');
end
