function tag = atlas_family_short_tag(name)
%ATLAS_FAMILY_SHORT_TAG  Compact tag for file names (Windows path-safe).

if isstring(name)
    name = char(name);
end
name = strtrim(name);

% Common cases
if contains(name, 'Lyapunov', 'IgnoreCase', true)
    if contains(name, 'L1', 'IgnoreCase', true)
        tag = 'LyapL1';
        return;
    elseif contains(name, 'L2', 'IgnoreCase', true)
        tag = 'LyapL2';
        return;
    elseif contains(name, 'L3', 'IgnoreCase', true)
        tag = 'LyapL3';
        return;
    end
end


% Fallback: sanitize + shorten
tag = regexprep(name, '[^A-Za-z0-9]+', '');
if isempty(tag)
    tag = 'fam';
end
if numel(tag) > 12
    tag = tag(1:12);
end
end
