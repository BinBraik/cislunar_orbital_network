function [minDV, dvlb, dvpatch, tof_at_minDV, voxelId_DV, ...
          minTOF, dv_at_minTOF, voxelId_TOF] = ...
        overlap_proxy_pair(FA, FB, grid3, cfg, VU_mps)
%OVERLAP_PROXY_PAIR  Compute DV-proxy and TOF-proxy winners for one pair,
% using pre-computed per-family voxel footprints (see overlap_proxy_footprint).
%
% For every candidate overlap voxel between FA and FB this scores two
% independent proxy costs:
%   dv_proxy  = dv_lb (sum of per-family min DV) + dv_patch (potential-based
%               plane-change estimate)
%   tof_proxy = t_min_A + t_min_B   (sum of per-family min TOF)
%
% and reports BOTH winners:
%   - the min-DV voxel   (minDV,  + TOF read off at that same voxel)
%   - the min-TOF voxel  (minTOF, + DV read off at that same voxel)
%
% This mirrors the original min-DV-only selection (which is preserved
% exactly as minDV/dvlb/dvpatch/tof_at_minDV/voxelId_DV) and adds an
% analogous, independently-optimized min-TOF selection.
%
% FA / FB are compact structs produced by overlap_proxy_footprint; they
% contain unique sorted voxel-ID vectors + per-voxel dv_min and t_min for
% FRS and BRS. No full row data is needed.
%
% minDV/dvlb/dvpatch/tof_at_minDV/voxelId_DV are numerically identical to
% the original local_run_pair's [minDV, dvlb, dvpatch, tof, voxelId].

minDV = NaN;  dvlb = NaN;  dvpatch = NaN;  tof_at_minDV = NaN;  voxelId_DV = NaN;
minTOF = NaN; dv_at_minTOF = NaN;          voxelId_TOF = NaN;
try
    % ── 1. Intersect FRS(A) with BRS(B) (both already sorted unique) ──────
    idsO = intersect(FA.uid_frs, FB.uid_brs);
    if isempty(idsO), return; end

    % ── 2. Unpack voxel grid indices ───────────────────────────────────────
    Ny = numel(grid3.y_centers);
    Nx = numel(grid3.x_centers);
    Nt = numel(grid3.th_centers);
    [iy, ix, ~] = ind2sub([Ny, Nx, Nt], idsO);

    % ── 3. Keep + primary-buffer filter (mirrors overlap_pair) ─────────
    bufFrac = 0.05;
    if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') ...
            && ~isempty(cfg.overlap.primary_buffer_frac)
        bufFrac = cfg.overlap.primary_buffer_frac;
    end
    if ~(isfield(cfg,'sys') && isfield(cfg.sys,'RE_nd') && isfield(cfg.sys,'RM_nd'))
        error('cfg.sys.RE_nd and cfg.sys.RM_nd required for primary buffer filter.');
    end
    RE = cfg.sys.RE_nd;
    RM = cfg.sys.RM_nd;
    mu = FA.mu;

    % Keep mask — all families share grid3_base so keepA = keepB = grid3.Keep
    if isfield(grid3,'Keep') && ~isempty(grid3.Keep)
        keepXY = logical(grid3.Keep);
        if ~isequal(size(keepXY), [Ny, Nx]), keepXY = keepXY.'; end
        okKeep = keepXY(sub2ind([Ny, Nx], iy, ix));
    else
        okKeep = true(numel(idsO), 1);
    end

    x = grid3.x_centers(ix);
    y = grid3.y_centers(iy);
    okEarth = hypot(x + mu, y)     > (1 + bufFrac) * RE;
    okMoon  = hypot(x - (1-mu), y) > (1 + bufFrac) * RM;
    ok = okKeep(:) & okEarth(:) & okMoon(:);

    idsO = idsO(ok);
    if isempty(idsO), return; end
    ix = ix(ok);  iy = iy(ok);

    % ── 4. Look up pre-computed per-voxel DV and TOF ──────────────────────
    % FA.uid_frs and FB.uid_brs are sorted → ismember is fast
    [~, locA] = ismember(idsO, FA.uid_frs);
    [~, locB] = ismember(idsO, FB.uid_brs);
    dv_min_A = FA.dv_min_frs(locA);
    dv_min_B = FB.dv_min_brs(locB);
    t_min_A  = FA.t_min_frs(locA);
    t_min_B  = FB.t_min_brs(locB);

    % ── 5. DV proxy (identical formula to original local_run_pair) ─────────
    x_ok = grid3.x_centers(ix);
    y_ok = grid3.y_centers(iy);
    CJstar = min(FA.CJ, FB.CJ);
    pot = cr3bp_potential(x_ok(:), y_ok(:), mu);
    v_box = sqrt(max(2 * pot.U - CJstar, 0));
    dv_patch_vec = 2 * v_box .* sin(abs(grid3.dtheta) / 2) * VU_mps;
    dv_lb_vec    = dv_min_A(:) + dv_min_B(:);
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    % ── 6. TOF proxy — min TOF per family, summed over the pair ────────────
    tof_proxy = t_min_A(:) + t_min_B(:);

    % ── 7. Min-DV winner (unchanged) ────────────────────────────────────────
    valid = isfinite(dv_proxy);
    if any(valid)
        idxValid  = find(valid);
        [~, iLoc] = min(dv_proxy(idxValid));
        iWin      = idxValid(iLoc);

        minDV        = dv_proxy(iWin);
        dvlb         = dv_lb_vec(iWin);
        dvpatch      = dv_patch_vec(iWin);
        voxelId_DV   = idsO(iWin);
        tof_at_minDV = tof_proxy(iWin);
    end

    % ── 8. Min-TOF winner — independent argmin over the same candidates ────
    validTOF = isfinite(tof_proxy);
    if any(validTOF)
        idxValidTOF   = find(validTOF);
        [~, iLocTOF]  = min(tof_proxy(idxValidTOF));
        iWinTOF       = idxValidTOF(iLocTOF);

        minTOF        = tof_proxy(iWinTOF);
        dv_at_minTOF  = dv_proxy(iWinTOF);
        voxelId_TOF   = idsO(iWinTOF);
    end
catch ME
    warning('[overlap_proxy_pair] %s', ME.message);
end
end
