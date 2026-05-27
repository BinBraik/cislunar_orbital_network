function [SeedsUpper, SeedsLower, meta] = atlas_seed_segments(xpo, ypo, thpo, y_eps, ds_seed, minSegPts)
% Supports multiple upper segments, robust interp1 via unique(s).

    meta = struct();

    isUpper = (ypo >= y_eps);

    segments = findContiguousSegments(isUpper);
    segments = mergeWrappedSegmentsIfNeeded(segments, isUpper);
    segments = segments(cellfun(@numel, segments) >= minSegPts);
    meta.nSegments = numel(segments);

    xU=[]; yU=[]; thU=[];

    for k=1:numel(segments)
        idx = segments{k};
        x_seg = xpo(idx); y_seg = ypo(idx); th_seg = thpo(idx);

        ds = hypot(diff(x_seg), diff(y_seg));
        s  = [0; cumsum(ds)];

        [sU, ia] = unique(s, 'stable');
        x_seg = x_seg(ia); y_seg = y_seg(ia); th_seg = th_seg(ia);

        if numel(sU) < 2, continue; end
        L = sU(end);
        if L <= 0, continue; end

        s_targets = 0:ds_seed:L;
        if s_targets(end) ~= L, s_targets = [s_targets, L]; end

        x_seed  = interp1(sU, x_seg,  s_targets, 'linear');
        y_seed  = interp1(sU, y_seg,  s_targets, 'linear');
        th_seed = interp1(sU, th_seg, s_targets, 'linear');

        xU  = [xU;  x_seed(:)]; %#ok<AGROW>
        yU  = [yU;  y_seed(:)]; %#ok<AGROW>
        thU = [thU; th_seed(:)]; %#ok<AGROW>
    end

    SeedsUpper = [xU(:), yU(:), thU(:)];
    SeedsLower = [xU(:), -yU(:), wrap_to_pi(pi - thU(:))];
end

function segments = findContiguousSegments(mask)
    idx = find(mask);
    segments = {};
    if isempty(idx), return; end
    breaks = [1; find(diff(idx) > 1) + 1; numel(idx) + 1];
    for k = 1:numel(breaks)-1
        segments{end+1} = idx(breaks(k):breaks(k+1)-1); %#ok<AGROW>
    end
end

function segments = mergeWrappedSegmentsIfNeeded(segments, mask)
    if isempty(segments), return; end
    if mask(1) && mask(end) && numel(segments) >= 2
        merged = [segments{end}; segments{1}];
        segments = [{merged}, segments(2:end-1)];
    end
end
