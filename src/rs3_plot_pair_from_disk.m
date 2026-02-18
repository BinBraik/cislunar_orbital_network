function rs3_plot_pair_from_disk(pairDir, varargin)
%RS3_PLOT_PAIR_FROM_DISK  Re-create Step 6/7/8 story figures from saved pair results.
%
% Usage
%   rs3_plot_pair_from_disk(pairDir)
%
% This is intended for Step 9 batch runs where plots are disabled.
% The function loads saved Step 6/7/8 MAT files from pairDir, reloads the
% cached family atlases, and regenerates the story figures (PNG + FIG).

if nargin < 1 || isempty(pairDir)
    error('pairDir is required.');
end
if ~exist(pairDir,'dir')
    error('pairDir not found: %s', pairDir);
end

% ---- Load cfg snapshot (pair-local preferred) ----
cfg = [];
pairCfgPath = fullfile(pairDir,'cfg_snapshot_pair.mat');
if exist(pairCfgPath,'file')
    S = load(pairCfgPath,'cfg_snapshot');
    cfg = S.cfg_snapshot;
else
    % try parent cfg_snapshot.mat (results/<tag>/cfg_snapshot.mat)
    parent = fileparts(fileparts(pairDir));
    parentCfgPath = fullfile(parent,'cfg_snapshot.mat');
    if exist(parentCfgPath,'file')
        S = load(parentCfgPath,'cfg_snapshot');
        cfg = S.cfg_snapshot;
    end
end
if isempty(cfg)
    cfg = rs3_cfg_defaults();
end

% enable plotting and .fig saving for inspection
cfg.io.save_figs = true;
cfg.io.fig_visible = 'on';
cfg.diag.story_plots = true;

% ---- Rebuild grid ----
grid3 = rs3_grid_make(cfg);
rs3_grid_validate(grid3, cfg);

% ---- Load step mats ----
step6file = local_find_one(pairDir, 'step6_overlap_*.mat');
S6 = load(step6file,'O');
O = S6.O;

step7file = local_find_zero_or_one(pairDir, 'step7_refine_*.mat');
R = [];
if ~isempty(step7file)
    S7 = load(step7file,'R');
    R = S7.R;
end

step8file = local_find_zero_or_one(pairDir, 'step8_score_*.mat');
P = [];
if ~isempty(step8file)
    S8 = load(step8file,'P');
    P = S8.P;
end

% ---- Load families from cache ----
Aname = O.A_name;
Bname = O.B_name;
A = rs3_prepare_or_load_family(Aname, cfg, grid3);
B = rs3_prepare_or_load_family(Bname, cfg, grid3);
KeepJoin = A.grid3.Keep & B.grid3.Keep;

% occupancy voxel ids for story plots
rowsA = A.Step4.rows;
rowsB = B.Step4.rows;
[vidA,~] = rs3_rows_to_voxel_ids(rowsA, grid3, KeepJoin);
[vidB,~] = rs3_rows_to_voxel_ids(rowsB, grid3, KeepJoin);

[poA, poAm] = rs3_get_po_xy(A, cfg);
[poB, poBm] = rs3_get_po_xy(B, cfg);

% ---- Step 6 story ----
if ~isempty(O) && isfield(O,'nOverlap') && O.nOverlap > 0
    rs3_step6_overlap_story(grid3, KeepJoin, vidA, vidB, O, cfg, pairDir, poA, poAm, poB, poBm);
end

% ---- Step 7 story ----
if ~isempty(R) && isfield(R,'Regions') && isfield(R,'stats')
    vidO = uint64(O.voxel_ids);
    rs3_step7_refine_story(grid3, KeepJoin, vidA, vidB, vidO, R.Regions, R.stats, cfg, pairDir, A.name, B.name, poA, poAm, poB, poBm);
end

% ---- Step 8 story + best trajectory ----
if ~isempty(P) && isfield(P,'best') && ~isempty(P.best)
    rs3_step8_score_story(P, grid3, KeepJoin, O, cfg, pairDir, poA, poAm, poB, poBm);
    rs3_step8_plot_best_trajectory(P, A, B, grid3, KeepJoin, O, cfg, pairDir, poA, poAm, poB, poBm);
end

fprintf('[rs3] Plots regenerated in: %s\n', pairDir);
end

function f = local_find_one(dir0, pat)
dd = dir(fullfile(dir0, pat));
if isempty(dd)
    error('File not found: %s', fullfile(dir0, pat));
end
f = fullfile(dd(1).folder, dd(1).name);
end

function f = local_find_zero_or_one(dir0, pat)
dd = dir(fullfile(dir0, pat));
if isempty(dd)
    f = '';
    return;
end
f = fullfile(dd(1).folder, dd(1).name);
end
