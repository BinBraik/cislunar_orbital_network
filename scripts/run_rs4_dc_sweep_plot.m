%% RUN_RS4_DC_SWEEP_PLOT
% Replay and save plots from a completed rs4_dc_sweep run.
%
% Reads rs4_dc_sweep_results.mat (produced by run_rs4_dc_sweep.m) and calls
% rs4_diffcorr_visualize for each selected pair, saving a side-by-side
% before/after PNG plus a plain-text info file per pair and a summary CSV.
%
% USER SETTINGS
%   RESULTS_FILE  : full path to rs4_dc_sweep_results.mat
%                   Leave '' to auto-detect newest run.
%
%   PLOT_PAIRS    : which pairs to plot —
%                     'all'              → every pair with valid trajectory data
%                     [1 3 7]            → specific pair indices (rows of pairs_ij)
%                     {'FamA','FamB'}    → single pair, looked up by name
%                     {{'FA1','FB1'};
%                      {'FA2','FB2'}}    → several named pairs
%
%   SKIP_DIVERGED : true  → skip pairs where DC did not converge
%                   false → plot all (converged status shown in title/info file)
%
% Outputs  (written to <run_root>/rs4_dc_sweep_plots/)
%   <pairSafe>/<pairSafe>_diffcorr_compare.png   — side-by-side before/after
%   <pairSafe>/<pairSafe>_info.txt               — plain-text metrics
%   transfer_summary.csv                          — one row per pair

clear; clc;

% ---- path setup ----------------------------------------------------------
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% =========================================================================
% USER SETTINGS
% =========================================================================
RESULTS_FILE  = '';      % leave '' to auto-detect
PLOT_PAIRS    = 'all';   % 'all' | [idx] | {'FamA','FamB'} | {{'FA1','FB1'}; ...}
SKIP_DIVERGED = false;   % true = skip unconverged pairs

% =========================================================================
% AUTO-DETECT results file
% =========================================================================
if isempty(RESULTS_FILE)
    candidates = dir(fullfile(repoRoot, 'rs3_results', '*', ...
                              'rs4_dc_sweep', 'rs4_dc_sweep_results.mat'));
    if isempty(candidates)
        error('[dc_plot] Cannot find rs4_dc_sweep_results.mat. Set RESULTS_FILE manually.');
    end
    [~, ord]     = sort([candidates.datenum], 'descend');
    RESULTS_FILE = fullfile(candidates(ord(1)).folder, candidates(ord(1)).name);
end
fprintf('[dc_plot] Results file : %s\n', RESULTS_FILE);

% =========================================================================
% LOAD RESULTS
% =========================================================================
res = load(RESULTS_FILE, ...
    'traj_cell', 'pairs_ij', 'families', ...
    'before_mat', 'after_mat', 'BEFORE_COLS', 'AFTER_COLS', 'cfg');

traj_cell   = res.traj_cell;
pairs_ij    = res.pairs_ij;
families    = res.families;
before_mat  = res.before_mat;
after_mat   = res.after_mat;
BEFORE_COLS = res.BEFORE_COLS;
AFTER_COLS  = res.AFTER_COLS;
cfg         = res.cfg;
nPairs      = size(pairs_ij, 1);
nFam        = numel(families);

% =========================================================================
% OUTPUT FOLDER  (sibling to rs4_dc_sweep/)
% =========================================================================
dcRoot    = fileparts(RESULTS_FILE);
plotsRoot = fullfile(fileparts(dcRoot), 'rs4_dc_sweep_plots');
if ~exist(plotsRoot, 'dir'), mkdir(plotsRoot); end
fprintf('[dc_plot] Output folder: %s\n\n', plotsRoot);

% =========================================================================
% RESOLVE WHICH PAIR INDICES TO PLOT
% =========================================================================
if ischar(PLOT_PAIRS) && strcmpi(PLOT_PAIRS, 'all')
    plot_idx = 1:nPairs;

elseif isnumeric(PLOT_PAIRS)
    plot_idx = PLOT_PAIRS(:)';

elseif iscell(PLOT_PAIRS)
    % {'FamA','FamB'}  or  {{'FA1','FB1'};{'FA2','FB2'}}
    if ischar(PLOT_PAIRS{1})
        named_pairs = {PLOT_PAIRS};   % wrap single pair
    else
        named_pairs = PLOT_PAIRS;
    end
    plot_idx = [];
    for k = 1:numel(named_pairs)
        fa = named_pairs{k}{1};
        fb = named_pairs{k}{2};
        ia = find(strcmp(families, fa), 1);
        ib = find(strcmp(families, fb), 1);
        if isempty(ia) || isempty(ib)
            warning('[dc_plot] Family not found: "%s" or "%s"', fa, fb);
            continue;
        end
        if ia > ib, [ia, ib] = deal(ib, ia); end   % enforce i < j
        row = find(pairs_ij(:,1)==ia & pairs_ij(:,2)==ib, 1);
        if isempty(row)
            warning('[dc_plot] Pair (%s, %s) not in results.', fa, fb);
            continue;
        end
        plot_idx(end+1) = row; %#ok<AGROW>
    end

else
    error('[dc_plot] Unrecognised PLOT_PAIRS format.');
end

% Drop pairs whose traj_cell entry is empty (extract or DC failed entirely)
has_traj = ~cellfun('isempty', traj_cell);
plot_idx = plot_idx(has_traj(plot_idx));

fprintf('[dc_plot] Pairs selected: %d\n', numel(plot_idx));
if numel(plot_idx) == 0
    fprintf('[dc_plot] Nothing to plot.\n');
    return;
end

% =========================================================================
% LOAD FAMILY STRUCTS + REBUILD Xpo
% =========================================================================
cfg.cache.rebuild = false;
grid3 = rs3_grid_make(cfg);
fprintf('[dc_plot] Loading %d family structs ...\n', nFam);

family_list = cell(nFam, 1);
for k = 1:nFam
    [family_list{k}, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    family_list{k} = local_ensure_xpo(family_list{k}, ...
        cfg.propag.relTol, cfg.propag.absTol, 1001);
end
fprintf('[dc_plot] Families loaded.\n\n');

% =========================================================================
% CONFIG FOR PLOTTING  (save every PNG, never show figures)
% =========================================================================
cfg_plot                 = cfg;
cfg_plot.io.save_figs    = true;
cfg_plot.io.save_fig     = false;   % no .fig files
cfg_plot.io.fig_visible  = 'off';

% =========================================================================
% MAIN PLOT LOOP
% =========================================================================
fprintf('[dc_plot] Plotting %d pairs ...\n\n', numel(plot_idx));

for kp = 1:numel(plot_idx)
    p  = plot_idx(kp);
    tr = traj_cell{p};

    ii = pairs_ij(p, 1);
    jj = pairs_ij(p, 2);

    if SKIP_DIVERGED && ~tr.converged
        fprintf('  [%3d/%d]  SKIP (diverged): %s → %s\n', ...
            kp, numel(plot_idx), families{ii}, families{jj});
        continue;
    end

    fprintf('  [%3d/%d]  %s → %s   conv=%d  DV: %.1f → %.1f m/s\n', ...
        kp, numel(plot_idx), families{ii}, families{jj}, ...
        tr.converged, tr.T_DV_total_mps, tr.DV_total_mps);

    SA = family_list{ii};
    SB = family_list{jj};

    % per-pair output subfolder
    pairTag  = sprintf('%s__TO__%s', tr.famA, tr.famB);
    pairSafe = rs3_sanitize_fname(pairTag);
    pairDir  = fullfile(plotsRoot, pairSafe);
    if ~exist(pairDir, 'dir'), mkdir(pairDir); end

    % Reconstruct T (before DC) and Tc (after DC) from traj struct
    T  = local_traj_to_T(tr);
    Tc = local_traj_to_Tc(tr);

    % Draw and save side-by-side figure
    try
        rs4_diffcorr_visualize(T, Tc, SA, SB, cfg_plot, pairDir, pairTag);
    catch ME
        warning('[dc_plot] Visualize failed for %s → %s: %s', ...
            families{ii}, families{jj}, ME.message);
    end
    close all;

    % Write plain-text info file
    local_write_info(tr, before_mat(p,:), after_mat(p,:), ...
        BEFORE_COLS, AFTER_COLS, pairDir, pairSafe);
end

% =========================================================================
% SUMMARY CSV  (all pairs, NaN where computation failed)
% =========================================================================
local_write_summary_csv(traj_cell, before_mat, after_mat, pairs_ij, ...
    families, BEFORE_COLS, AFTER_COLS, plotsRoot);

fprintf('\n[dc_plot] Done.  Plots saved in:\n  %s\n', plotsRoot);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function T = local_traj_to_T(tr)
% Reconstruct the before-DC T struct from a traj_cell entry.
T = struct();
T.IC_A              = tr.T_IC_A;
T.IC_B_frs          = tr.T_IC_B_frs;
T.seed_A            = tr.T_seed_A;
T.seed_B_frs        = tr.T_seed_B_frs;
T.XA                = tr.T_XA;
T.x_B               = tr.T_x_B;
T.y_B               = tr.T_y_B;
T.th_B              = tr.T_th_B;
T.DV_turn_A_mps     = tr.T_DV_turn_A_mps;
T.DV_patch_mps      = tr.T_DV_patch_mps;
T.DV_turn_B_mps     = tr.T_DV_turn_B_mps;
T.DV_total_true_mps = tr.T_DV_total_mps;
T.tof_A_days        = tr.T_tof_A_days;
T.tof_B_days        = tr.T_tof_B_days;
T.vid               = tr.vid;
T.xc                = tr.T_xc;
T.yc                = tr.T_yc;
T.i_star            = tr.T_i_star;
T.j_star            = tr.T_j_star;
T.dA_nd             = tr.T_dA_nd;
T.dB_nd             = tr.T_dB_nd;
T.delta_th_rad      = tr.T_delta_th_rad;
T.DV_proxy_mps      = tr.T_DV_proxy_mps;
end

% -------------------------------------------------------------------------

function Tc = local_traj_to_Tc(tr)
% Reconstruct the after-DC Tc struct from a traj_cell entry.
Tc = struct();
Tc.IC_A          = tr.IC_A;
Tc.IC_B_frs      = tr.IC_B_frs;
Tc.seed_A        = tr.seed_A;
Tc.seed_B_frs    = tr.seed_B_frs;
Tc.XA            = tr.XA;
Tc.x_B           = tr.x_B;
Tc.y_B           = tr.y_B;
Tc.th_B          = tr.th_B;
Tc.xp            = tr.xp;
Tc.yp            = tr.yp;
Tc.DV_turn_A_mps = tr.DV_turn_A_mps;
Tc.DV_patch_mps  = tr.DV_patch_mps;
Tc.DV_turn_B_mps = tr.DV_turn_B_mps;
Tc.DV_total_mps  = tr.DV_total_mps;
Tc.tof_A_days    = tr.tof_A_days;
Tc.tof_B_days    = tr.tof_B_days;
Tc.delta_th_rad  = tr.delta_th_rad;
Tc.r_final       = tr.r_final;
Tc.exitflag      = tr.exitflag;
Tc.converged     = tr.converged;
end

% -------------------------------------------------------------------------

function local_write_info(tr, bef_row, aft_row, BEFORE_COLS, AFTER_COLS, pairDir, pairSafe)
% Write a plain-text metrics file alongside the PNG.
fid = fopen(fullfile(pairDir, [pairSafe '_info.txt']), 'w');
if fid < 0, return; end

gb = @(name) bef_row(find(strcmp(BEFORE_COLS, name), 1));
ga = @(name) aft_row(find(strcmp(AFTER_COLS,  name), 1));

fprintf(fid, '=== Transfer: %s  -->  %s ===\n\n', tr.famA, tr.famB);
fprintf(fid, 'Voxel ID       : %d\n', tr.vid);

fprintf(fid, '\n--- Before DC ---\n');
fprintf(fid, 'DV_turn_A      : %8.2f  m/s\n',  gb('DV_turn_A_mps'));
fprintf(fid, 'DV_patch       : %8.2f  m/s\n',  gb('DV_patch_mps'));
fprintf(fid, 'DV_turn_B      : %8.2f  m/s\n',  gb('DV_turn_B_mps'));
fprintf(fid, 'DV_total       : %8.2f  m/s\n',  gb('DV_total_mps'));
fprintf(fid, 'TOF_A          : %8.3f  days\n', gb('tof_A_days'));
fprintf(fid, 'TOF_B          : %8.3f  days\n', gb('tof_B_days'));
fprintf(fid, 'delta_theta    : %8.4f  deg\n',  gb('delta_th_deg'));
fprintf(fid, 'miss_A         : %8.6f  nd\n',   gb('dA_nd'));
fprintf(fid, 'miss_B         : %8.6f  nd\n',   gb('dB_nd'));

fprintf(fid, '\n--- After DC ---\n');
fprintf(fid, 'DV_turn_A      : %8.2f  m/s\n',  ga('DV_turn_A_mps'));
fprintf(fid, 'DV_patch       : %8.2f  m/s\n',  ga('DV_patch_mps'));
fprintf(fid, 'DV_turn_B      : %8.2f  m/s\n',  ga('DV_turn_B_mps'));
fprintf(fid, 'DV_total       : %8.2f  m/s\n',  ga('DV_total_mps'));
fprintf(fid, 'TOF_A          : %8.3f  days\n', ga('tof_A_days'));
fprintf(fid, 'TOF_B          : %8.3f  days\n', ga('tof_B_days'));
fprintf(fid, 'delta_theta    : %8.4f  deg\n',  ga('delta_th_deg'));
fprintf(fid, 'exitflag       : %d\n',           ga('exitflag'));
fprintf(fid, 'converged      : %d\n',           ga('converged'));
fprintf(fid, 'r_norm_scaled  : %.4e\n',         ga('r_norm_scaled'));

dv_improve = gb('DV_total_mps') - ga('DV_total_mps');
fprintf(fid, '\nDV improvement : %+.2f  m/s\n', dv_improve);

fclose(fid);
end

% -------------------------------------------------------------------------

function local_write_summary_csv(traj_cell, before_mat, after_mat, ...
        pairs_ij, families, BEFORE_COLS, AFTER_COLS, plotsRoot)
% Write one CSV row per pair (NaN where computation failed).
csvFile = fullfile(plotsRoot, 'transfer_summary.csv');
fid = fopen(csvFile, 'w');
if fid < 0
    warning('[dc_plot] Cannot write summary CSV: %s', csvFile);
    return;
end

% Header
fprintf(fid, ['pair_idx,fam_A,fam_B,' ...
    'bef_DV_turn_A,bef_DV_patch,bef_DV_turn_B,bef_DV_total,' ...
    'bef_tof_A_days,bef_tof_B_days,bef_delta_th_deg,bef_dA_nd,bef_dB_nd,' ...
    'aft_DV_turn_A,aft_DV_patch,aft_DV_turn_B,aft_DV_total,' ...
    'aft_tof_A_days,aft_tof_B_days,aft_delta_th_deg,' ...
    'exitflag,converged,r_norm_scaled,DV_improvement\n']);

bc = @(name) find(strcmp(BEFORE_COLS, name), 1);
ac = @(name) find(strcmp(AFTER_COLS,  name), 1);

for p = 1:size(pairs_ij, 1)
    ii = pairs_ij(p, 1);
    jj = pairs_ij(p, 2);
    b  = before_mat(p, :);
    a  = after_mat(p,  :);

    dv_improve = b(bc('DV_total_mps')) - a(ac('DV_total_mps'));

    fprintf(fid, '%d,%s,%s,', p, families{ii}, families{jj});
    fprintf(fid, '%.3f,%.3f,%.3f,%.3f,', ...
        b(bc('DV_turn_A_mps')), b(bc('DV_patch_mps')), ...
        b(bc('DV_turn_B_mps')), b(bc('DV_total_mps')));
    fprintf(fid, '%.4f,%.4f,%.4f,%.6f,%.6f,', ...
        b(bc('tof_A_days')), b(bc('tof_B_days')), ...
        b(bc('delta_th_deg')), b(bc('dA_nd')), b(bc('dB_nd')));
    fprintf(fid, '%.3f,%.3f,%.3f,%.3f,', ...
        a(ac('DV_turn_A_mps')), a(ac('DV_patch_mps')), ...
        a(ac('DV_turn_B_mps')), a(ac('DV_total_mps')));
    fprintf(fid, '%.4f,%.4f,%.4f,', ...
        a(ac('tof_A_days')), a(ac('tof_B_days')), a(ac('delta_th_deg')));
    fprintf(fid, '%d,%d,%.4e,%.3f\n', ...
        a(ac('exitflag')), a(ac('converged')), a(ac('r_norm_scaled')), dv_improve);
end

fclose(fid);
fprintf('[dc_plot] Summary CSV  → %s\n', csvFile);
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
% Re-integrate the periodic orbit if Xpo was stripped from cache.
    if isfield(S,'Xpo') && ~isempty(S.Xpo) && ...
       isfield(S,'t_dense') && ~isempty(S.t_dense)
        return;
    end
    fprintf('    [ensure_xpo] rebuilding Xpo for "%s" ...\n', S.name);
    opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
    solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                     [0, S.Tf_PO], S.X0, opts);
    t_dense = linspace(0, S.Tf_PO, N_po)';
    S.t_dense = t_dense;
    S.Xpo     = deval(solPO, t_dense)';
end
