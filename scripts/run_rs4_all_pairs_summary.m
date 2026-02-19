%% RUN_RS4_ALL_PAIRS_SUMMARY
% Batch overlap runner over all family pairs with per-pair top-10 voxel export
% and summary CSV outputs.
clear; clc;

% --- hard-pin this repo paths ---
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot,'src'));
addpath(fullfile(repoRoot,'scripts'));
rehash;

cfg = rs3_cfg_defaults();

% ===================== USER KNOBS =====================
% If you want a specific subset/order, edit this list.
families = { ...
    'Lyapunov L1', ...
    'Lyapunov L2', ...
    'Cycler 21', ...
    'Cycler 11a', ...
    'Cycler 11b', ...
    'Cycler 32', ...
    'Resonant 2to1 Stable', ...
    'Resonant 2to1 Unstable', ...
    'Resonant 3to1 Stable', ...
    'Resonant 3to1 Unstable', ...
    'Resonant 5to2 Stable', ...
    'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit' ...
    };

cfg.families.list = families;
cfg.families.test_only = false;

cfg.io.save_figs   = true;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% grid settings (must match for both atlases)
cfg.grid.dx     = 0.01;
cfg.grid.dy     = 0.01;
cfg.grid.dtheta = deg2rad(2);

cfg.seed.ds_seed   = 0.02;

% propagation/fan
cfg.propag.Tmax    = pi;
cfg.fan.DV_cap_nd  = 0.2;
cfg.fan.dtheta_fan = deg2rad(1.0);
cfg.propag.absTol  = 1e-9;
cfg.propag.relTol  = 1e-9;
cfg.propag.v2tol   = 1e-8;
cfg.log.step_len_factor = 0.5;
cfg.log.maxstep_factor  = 0.5;

% zoom optional
cfg.diag.zoom.enable = false;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.45 0.45];

% Figure toggles (optional)
cfg.plot.rs4.overlap_xy   = true;
cfg.plot.rs4.overlap_xyz  = false;
cfg.plot.rs4.combo_xy     = true;
cfg.plot.rs4.combo_xyz    = false;
cfg.plot.rs4.bounds_lb    = true;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = true;

rs3_cfg_validate(cfg);

batchTag = sprintf('rs4_pairs_%dfam', numel(families));
outRoot = fullfile(cfg.io.out_root, cfg.io.tag, batchTag);
if ~exist(outRoot,'dir'), mkdir(outRoot); end
summaryDir = fullfile(outRoot, 'Summary');
if ~exist(summaryDir,'dir'), mkdir(summaryDir); end

% ===================== BUILD GRID =====================
grid3 = rs3_grid_make(cfg);
if exist('rs3_grid_validate','file')==2
    rs3_grid_validate(grid3, cfg);
end

% ===================== PREBUILD / LOAD ALL FAMILIES =====================
N = numel(families);
Sall = cell(N,1);
InfoAll = cell(N,1);
for i = 1:N
    fam = families{i};
    fprintf('[rs4-batch] build/load family %d/%d: %s\n', i, N, fam);
    [Sall{i}, InfoAll{i}] = rs3_prepare_or_load_family(fam, cfg, grid3);
end

% ===================== SUMMARY HOLDERS =====================
minDVproxyMat = nan(N,N);

nPairs = N*(N-1)/2;
colFamA = cell(nPairs,1);
colFamB = cell(nPairs,1);
colMinDV = nan(nPairs,1);
colDVlb = nan(nPairs,1);
colDVpatch = nan(nPairs,1);
colTOF = nan(nPairs,1);
colVoxelId = nan(nPairs,1);

pairRow = 0;

% ===================== PAIR LOOP =====================
for i = 1:N
    for j = i+1:N
        pairRow = pairRow + 1;
        famA = families{i};
        famB = families{j};

        colFamA{pairRow} = famA;
        colFamB{pairRow} = famB;

        pairTag = sprintf('%s__TO__%s', famA, famB);
        pairSafe = rs3_sanitize_fname(pairTag);
        pairDir = fullfile(outRoot, pairSafe);
        if ~exist(pairDir,'dir'), mkdir(pairDir); end

        fprintf('\n[rs4-batch] pair %d/%d: %s -> %s\n', pairRow, nPairs, famA, famB);

        SA = Sall{i};
        SB = Sall{j};

        try
            O = rs4_overlap_pair(SA, SB, cfg);
            save(fullfile(pairDir, ['rs4_' pairSafe '_overlap.mat']), 'O', '-v7.3');

            if isempty(O.ids)
                fprintf('[rs4-batch] no overlap for %s | %s\n', famA, famB);
                status = struct('has_overlap', false, 'familyA', famA, 'familyB', famB, ...
                    'min_dvproxy_mps', NaN, 'dv_lb_mps', NaN, 'dv_patch_ub_mps', NaN, ...
                    'tof_est_days', NaN, 'voxel_id', NaN);
                save(fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']), 'status', '-v7.3');
                continue;
            end

            V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);

            rs4_overlap_visualize(O, SA, SB, cfg, pairDir, pairTag);
            rs4_overlap_visualize_combo(SA, SB, O, cfg, pairDir, pairTag);
            B = rs4_overlap_visualize_bounds(V, SA, SB, cfg, pairDir, pairTag);

            % Winner selection by DVproxy
            dv = B.dv_proxy(:);
            idxFinite = find(isfinite(dv));

            if isempty(idxFinite)
                fprintf('[rs4-batch] no finite DVproxy for %s | %s\n', famA, famB);
                status = struct('has_overlap', true, 'familyA', famA, 'familyB', famB, ...
                    'min_dvproxy_mps', NaN, 'dv_lb_mps', NaN, 'dv_patch_ub_mps', NaN, ...
                    'tof_est_days', NaN, 'voxel_id', NaN);
                save(fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']), 'status', 'B', '-v7.3');
                continue;
            end

            [~, orderLocal] = sort(dv(idxFinite), 'ascend');
            iWinner = idxFinite(orderLocal(1));
            minDV = dv(iWinner);
            dvlbWinner = local_safe_index(B.dv_lb, iWinner);
            dvpatchWinner = local_safe_index(B.dv_patch_ub, iWinner);

            % Winner voxel info from flat V arrays
            tofWinner = NaN;
            voxelIdWinner = NaN;
            if iWinner >= 1 && iWinner <= numel(V.ids)
                voxelIdWinner = V.ids(iWinner);
                tofWinner = V.t_days_mean_A(iWinner) + V.t_days_mean_B(iWinner);
            end

            % Fill summary matrix and row
            minDVproxyMat(i,j) = minDV;
            minDVproxyMat(j,i) = minDV;

            colMinDV(pairRow) = minDV;
            colDVlb(pairRow) = dvlbWinner;
            colDVpatch(pairRow) = dvpatchWinner;
            colTOF(pairRow) = tofWinner;
            colVoxelId(pairRow) = voxelIdWinner;

            % Save winner + full bounds (DVproxy for all voxels) per pair
            winnerMeta = struct();
            winnerMeta.familyA = famA;
            winnerMeta.familyB = famB;
            winnerMeta.pairTag = pairTag;
            winnerMeta.generated = datestr(now, 31);
            winnerMeta.units = struct('dv','m/s','tof','days');
            winnerMeta.grid = struct('dx',grid3.dx,'dy',grid3.dy,'dtheta',grid3.dtheta);
            winnerMeta.iWinner = iWinner;
            winnerMeta.voxel_id = voxelIdWinner;
            winnerMeta.min_dvproxy_mps = minDV;
            winnerMeta.dv_lb_mps = dvlbWinner;
            winnerMeta.dv_patch_ub_mps = dvpatchWinner;
            winnerMeta.tof_est_days = tofWinner;

            save(fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']), ...
                'winnerMeta', 'B', '-v7.3');

            fprintf('[rs4-batch] min DVproxy = %.3f m/s (voxel id %d)\n', minDV, round(voxelIdWinner));

        catch ME
            warning('[rs4-batch] pair failed for %s | %s: %s', famA, famB, ME.message);
            status = struct('has_overlap', false, 'familyA', famA, 'familyB', famB, ...
                'error', ME.message, 'min_dvproxy_mps', NaN, 'dv_lb_mps', NaN, ...
                'dv_patch_ub_mps', NaN, 'tof_est_days', NaN, 'voxel_id', NaN);
            save(fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']), 'status', '-v7.3');
        end
    end
end

% ===================== SUMMARY CSV #1: matrix =====================
matrixCell = cell(N+1, N+1);
matrixCell{1,1} = 'Family';
for i = 1:N
    matrixCell{1, i+1} = families{i};
    matrixCell{i+1, 1} = families{i};
    for j = 1:N
        matrixCell{i+1, j+1} = minDVproxyMat(i,j);
    end
end
writecell(matrixCell, fullfile(summaryDir, 'minDVproxy_matrix.csv'));

% ===================== SUMMARY CSV #2: long table (top-1/pair) =====================
T = table(colFamA, colFamB, colMinDV, colDVlb, colDVpatch, colTOF, colVoxelId, ...
    'VariableNames', {'FamilyA','FamilyB','minDVproxy_mps','DVlb_mps','DVpatch_ub_mps','EstimatedTOF_days','VoxelId'});
writetable(T, fullfile(summaryDir, 'pair_winners_top1.csv'));

save(fullfile(summaryDir, 'batch_summary_workspace.mat'), 'families', 'minDVproxyMat', 'T', 'cfg', 'InfoAll', '-v7.3');

fprintf('\n[rs4-batch] done.\n');
fprintf('  root:    %s\n', outRoot);
fprintf('  summary: %s\n', summaryDir);

% ===================== local helpers =====================
function v = local_safe_index(arr, idx)
v = NaN;
if isempty(arr), return; end
arr = arr(:);
if idx >= 1 && idx <= numel(arr)
    v = arr(idx);
end
end

function v = local_cfg_get(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k), return; end
        cur = cur.(k);
    end
    if ~isempty(cur), v = cur; end
catch
    v = defaultVal;
end
end

