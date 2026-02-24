function st = rs3_cache_family_stats(S, cfg)
%RS3_CACHE_FAMILY_STATS  Compute basic stats for cache validation and size estimates.

st = struct();

% seeds/headings
st.Nseeds = size(S.SeedsUpper,1);
if isfield(S,'Step4') && isfield(S.Step4,'delta_lists')
    nh = cellfun(@numel, S.Step4.delta_lists);
    st.Nheads_total = sum(nh);
    st.Nheads_max = max(nh);
    st.Nheads_mean = mean(nh);
else
    st.Nheads_total = NaN;
    st.Nheads_max = NaN;
    st.Nheads_mean = NaN;
end

% hit table sizes
if isfield(S,'Step4')
    st.rows_FRS_upper = rs3_rows_count(S.Step4.rows_FRS_upper);
    st.rows_FRS_lower = rs3_rows_count(S.Step4.rows_FRS_lower);
    % BRS is R(FRS): BRS_upper = R(FRS_lower), BRS_lower = R(FRS_upper)
    st.rows_BRS_upper = st.rows_FRS_lower;  % same count (kept for cache inspect compat)
    st.rows_BRS_lower = st.rows_FRS_upper;

    st.uv_FRS_upper = numel(unique(rs3_rows_to_voxelid(S.Step4.rows_FRS_upper, S.grid3)));
    st.uv_FRS_lower = numel(unique(rs3_rows_to_voxelid(S.Step4.rows_FRS_lower, S.grid3)));
    st.uv_BRS_upper = st.uv_FRS_lower;     % kept for cache inspect compat
else
    st.rows_FRS_upper = NaN;
    st.rows_FRS_lower = NaN;
    st.rows_BRS_upper = NaN;
    st.rows_BRS_lower = NaN;
    st.uv_FRS_upper = NaN;
    st.uv_FRS_lower = NaN;
    st.uv_BRS_upper = NaN;
end

% rough size estimate (bytes)
if isfield(S,'Step4')
    nrows = st.rows_FRS_upper + st.rows_FRS_lower;
    % Packed format: ~14 bytes/row (not 64)
    isPacked = isfield(S.Step4,'packed') && S.Step4.packed;
    bpr = 14; if ~isPacked, bpr = 64; end
    st.bytes_est = nrows * bpr;
else
    st.bytes_est = 0;
end

st.cfg_version_tag = cfg.cache.version_tag;
end
