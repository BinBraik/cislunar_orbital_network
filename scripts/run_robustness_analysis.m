%% RUN_ROBUSTNESS_ANALYSIS
% Cross-run robustness analysis for the verification grid sweep.
% Loads minDVproxyMat from each run in rs4_verification_runs/ and computes:
%
%   Table A  — per-run graph summary
%              (edges, density, LCC size, avg shortest-path DV, reciprocity)
%
%   Table B  — per-family centrality per run
%              (degree, strength, harmonic closeness, betweenness, PageRank)
%
%   Table C  — cross-run robustness
%              (mean rank, std rank, best rank, worst rank, top-3 / top-5 count)
%              separately for: degree, strength, harmonic closeness, betweenness, PageRank
%
%   Table D  — edge persistence
%              (for each family pair: n_runs with edge, mean DV, std DV)
%
%   Figure 1 — Rank heatmap (families × runs, one subplot per metric)
%   Figure 2 — Top-5 families per metric with std-rank error bars
%   Figure 3 — Edge-persistence heatmap (N×N, colour = fraction of runs)
%   Figure 4 — Baseline network (Run 0)
%   Figure 5 — Coarse vs fine network comparison
%
% Outputs written to rs4_verification_runs/robustness_analysis/

clear; clc;

% ── repo paths ───────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ════════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ════════════════════════════════════════════════════════════════════════════

% Run IDs are auto-discovered: every subfolder of rs4_verification_runs/ that
% contains a result.mat is loaded.  Test runs added via run_verification_grid_sweep.m
% are picked up automatically — no manual editing needed here.
verifyRoot_probe = fullfile(repoRoot, 'rs4_verification_runs');
d = dir(fullfile(verifyRoot_probe, '*/result.mat'));
if isempty(d)
    error('[robustness] No result.mat files found under %s\nRun run_verification_grid_sweep.m first.', ...
        verifyRoot_probe);
end
run_ids = arrayfun(@(x) x.folder(length(verifyRoot_probe)+2:end), d, 'UniformOutput', false);

% Network parameters (must match what was used in the sweep runs)
DV_CAP_ND     = 0.2;   % non-dimensional DV cap
TMAX_ND       = pi;    % non-dimensional Tmax
BUDGET_FACTOR = 2;     % departure + arrival (same as run_network_centrality_sweep.m)

% Edge existence threshold: edge (i,j) kept if DV(i,j) ≤ DV_budget
% DV_budget = DV_CAP_ND * VU_mps * BUDGET_FACTOR  (computed below per run)

SAVE_FIGS = true;   % true = save PNG + FIG; false = display only

% ── Paths ────────────────────────────────────────────────────────────────────
verifyRoot = fullfile(repoRoot, 'rs4_verification_runs');
outDir     = fullfile(verifyRoot, 'robustness_analysis');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% ── Family names ─────────────────────────────────────────────────────────────
families_default = { ...
    'Lyapunov L1', 'Lyapunov L2', ...
    'Cycler 21', 'Cycler 11a', 'Cycler 11b', 'Cycler 32', ...
    'Resonant 2to1 Stable', 'Resonant 2to1 Unstable', ...
    'Resonant 3to1 Stable', 'Resonant 3to1 Unstable', ...
    'Resonant 5to2 Stable', 'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit' ...
    };

% ── Metric metadata ──────────────────────────────────────────────────────────
metric_ids     = {'degree', 'strength', 'harmonic_closeness', 'betweenness', 'pagerank'};
metric_labels  = {'Degree', 'Strength (m/s)', 'Harmonic Closeness', 'Betweenness', 'PageRank'};
nMetrics       = numel(metric_ids);

% ════════════════════════════════════════════════════════════════════════════
%  LOAD RESULTS
% ════════════════════════════════════════════════════════════════════════════
nRuns = numel(run_ids);
DV_all  = cell(nRuns, 1);
TOF_all = cell(nRuns, 1);
cfg_all = cell(nRuns, 1);
valid   = false(nRuns, 1);
families = families_default;
N = numel(families);

fprintf('[robustness] Loading %d run result files...\n', nRuns);
for r = 1:nRuns
    fpath = fullfile(verifyRoot, run_ids{r}, 'result.mat');
    if ~exist(fpath, 'file')
        warning('[robustness] Missing: %s  — skipping run.', fpath);
        continue;
    end
    d = load(fpath, 'minDVproxyMat', 'families', 'cfg');
    if isfield(d, 'TOFmat')
        d2 = load(fpath, 'TOFmat');
        TOF_all{r} = d2.TOFmat;
    end
    DV_all{r}  = d.minDVproxyMat;
    cfg_all{r} = d.cfg;
    % Use family list from first valid run
    if isfield(d,'families') && ~isempty(d.families)
        families = d.families;
        N = numel(families);
    end
    valid(r) = true;
    fprintf('[robustness]   loaded %s\n', run_ids{r});
end

valid_idx = find(valid);
nValid    = numel(valid_idx);
fprintf('[robustness] %d/%d runs loaded successfully.\n\n', nValid, nRuns);

if nValid == 0
    error('[robustness] No runs found. Run run_verification_grid_sweep.m first.');
end

% Short family names for figures
if exist('net_family_short_names', 'file') == 2
    short_names = net_family_short_names(families);
else
    short_names = cellfun(@(s) strtrim(s(end-min(7,numel(s)-1):end)), ...
        families, 'UniformOutput', false);
end

% ════════════════════════════════════════════════════════════════════════════
%  PRE-ALLOCATE STORAGE
% ════════════════════════════════════════════════════════════════════════════
deg_all = nan(N, nRuns);
str_all = nan(N, nRuns);
hc_all  = nan(N, nRuns);
btw_all = nan(N, nRuns);
pgr_all = nan(N, nRuns);

n_edges_all  = nan(1, nRuns);
density_all  = nan(1, nRuns);
lcc_all      = nan(1, nRuns);
avg_path_all = nan(1, nRuns);

% Edge-persistence accumulator [N×N upper-triangle counters, sum DV, sum DV²]
edge_count  = zeros(N, N);
edge_sumDV  = zeros(N, N);
edge_sumDV2 = zeros(N, N);

% ════════════════════════════════════════════════════════════════════════════
%  PER-RUN NETWORK METRICS
% ════════════════════════════════════════════════════════════════════════════
cfg_base = rs3_cfg_defaults();
VU_mps   = cfg_base.units.VU_mps;
TU_days  = cfg_base.units.TU_days;

DV_budget = DV_CAP_ND * VU_mps * BUDGET_FACTOR;   % m/s threshold for edge existence

fprintf('[robustness] Computing per-run metrics (DV budget = %.0f m/s)...\n', DV_budget);

for r = valid_idx'
    DV = DV_all{r};
    if isempty(DV), continue; end

    % Adjacency: edge exists if finite and within budget
    A = isfinite(DV) & (DV <= DV_budget);
    A = A & ~eye(N);   % no self-loops

    edges_kept = sum(sum(triu(A, 1)));
    if edges_kept == 0
        fprintf('[robustness]   Run %-22s: 0 edges — skipping metrics.\n', run_ids{r});
        continue;
    end

    % Weighted adjacency (0 for missing/infeasible, DV for feasible edges)
    W = DV .* A;   % 0 where no edge

    % Floyd-Warshall — W needs Inf for missing edges (not 0)
    W_fw = DV;
    W_fw(~A) = Inf;
    W_fw(1:N+1:end) = 0;   % zero diagonal
    if exist('net_floyd_warshall', 'file') == 2
        dist = net_floyd_warshall(W_fw);
    else
        dist = local_floyd_warshall(W_fw, N);
    end

    % LCC size
    if exist('net_lcc', 'file') == 2
        [~, lcc_sz] = net_lcc(A);
    else
        lcc_sz = local_lcc(A, N);
    end

    % ── Graph-level metrics ────────────────────────────────────────────────
    n_edges_all(r)  = edges_kept;
    density_all(r)  = edges_kept / (N*(N-1)/2);
    lcc_all(r)      = lcc_sz;
    off_d           = dist(~eye(N,'logical'));
    avg_path_all(r) = mean(off_d(isfinite(off_d)));

    % ── Degree ─────────────────────────────────────────────────────────────
    deg_all(:, r) = sum(A, 2);   % undirected: sum of row = number of neighbours

    % ── Strength (sum of edge DV costs incident to each node) ──────────────
    str_all(:, r) = sum(W, 2);

    % ── Harmonic closeness ─────────────────────────────────────────────────
    % hc(i) = (1/(N-1)) * Σ_{j≠i, finite} 1/dist(i,j)
    inv_dist = 1 ./ dist;
    inv_dist(isinf(dist) | isnan(dist)) = 0;
    inv_dist(1:N+1:end) = 0;   % exclude self
    hc_all(:, r) = sum(inv_dist, 2) / (N-1);

    % ── Betweenness (weighted, using MATLAB graph object) ──────────────────
    [ii_e, jj_e] = find(triu(A, 1));
    dv_e = DV(sub2ind([N,N], ii_e, jj_e));
    try
        G   = graph(ii_e, jj_e, dv_e, N);
        btw = centrality(G, 'betweenness', 'Cost', G.Edges.Weight);
    catch
        btw = local_betweenness(A, dist, N);
    end
    mx = max(btw);
    btw_all(:, r) = btw / max(mx, eps);   % normalise to [0,1]

    % ── PageRank ───────────────────────────────────────────────────────────
    try
        % Importance weight = 1/DV (stronger connection = easier to follow)
        imp = 1 ./ max(dv_e, 1);   % avoid /0; DV always > 0 for feasible edges
        G2  = graph(ii_e, jj_e, imp, N);
        pgr = centrality(G2, 'pagerank', 'FollowProbability', 0.85, ...
                         'Importance', G2.Edges.Weight);
    catch
        pgr = nan(N, 1);
    end
    pgr_all(:, r) = pgr;

    % ── Edge-persistence accumulator ───────────────────────────────────────
    for p = 1:numel(ii_e)
        i = ii_e(p);  j = jj_e(p);
        edge_count(i,j)  = edge_count(i,j)  + 1;
        edge_sumDV(i,j)  = edge_sumDV(i,j)  + dv_e(p);
        edge_sumDV2(i,j) = edge_sumDV2(i,j) + dv_e(p)^2;
    end

    fprintf('[robustness]   Run %-22s: %d edges  LCC=%d  avg_path=%.0f m/s\n', ...
        run_ids{r}, edges_kept, lcc_sz, avg_path_all(r));
end

% ════════════════════════════════════════════════════════════════════════════
%  TABLE A — Per-run graph summary
% ════════════════════════════════════════════════════════════════════════════
fprintf('\n[robustness] Building Table A (graph summary)...\n');
run_labels = run_ids(:);
T_A = table(run_labels, n_edges_all(:), density_all(:), lcc_all(:), avg_path_all(:), ...
    'VariableNames', {'RunID','Edges','Density','LCC_size','AvgPathDV_mps'});
writetable(T_A, fullfile(outDir, 'TableA_graph_summary.csv'));
disp(T_A);

% ════════════════════════════════════════════════════════════════════════════
%  TABLE B — Per-family centrality per run (one CSV per metric)
% ════════════════════════════════════════════════════════════════════════════
fprintf('[robustness] Building Table B (per-family centrality)...\n');
metric_data = {deg_all, str_all, hc_all, btw_all, pgr_all};

for m = 1:nMetrics
    dat = metric_data{m};           % [N × nRuns]
    col_names = ['Family', run_ids(:)'];
    T_B = array2table([families(:), num2cell(dat)], 'VariableNames', col_names);
    fname = sprintf('TableB_%s_per_run.csv', metric_ids{m});
    writetable(T_B, fullfile(outDir, fname));
end

% ════════════════════════════════════════════════════════════════════════════
%  RANK ARRAYS  (rank 1 = best / highest metric value)
% ════════════════════════════════════════════════════════════════════════════
rank_all = cell(nMetrics, 1);
for m = 1:nMetrics
    dat = metric_data{m};
    rnk = nan(N, nRuns);
    for r = valid_idx'
        col = dat(:, r);
        if all(isnan(col)), continue; end
        [~, order] = sort(col, 'descend', 'MissingPlacement', 'last');
        rnk(order, r) = 1:N;
    end
    rank_all{m} = rnk;
end

% ════════════════════════════════════════════════════════════════════════════
%  TABLE C — Cross-run robustness (one CSV per metric)
% ════════════════════════════════════════════════════════════════════════════
fprintf('[robustness] Building Table C (cross-run robustness)...\n');
for m = 1:nMetrics
    rnk    = rank_all{m}(:, valid_idx);   % only valid runs
    m_rank = nanmean(rnk, 2);
    s_rank = nanstd(rnk, 0, 2);
    b_rank = nanmin(rnk, [], 2);
    w_rank = nanmax(rnk, [], 2);
    top3   = sum(rnk <= 3, 2);
    top5   = sum(rnk <= 5, 2);

    T_C = table(families(:), m_rank, s_rank, b_rank, w_rank, top3, top5, ...
        'VariableNames', {'Family','MeanRank','StdRank','BestRank','WorstRank', ...
                          'Top3Count','Top5Count'});
    T_C = sortrows(T_C, 'MeanRank');

    fname = sprintf('TableC_%s_robustness.csv', metric_ids{m});
    writetable(T_C, fullfile(outDir, fname));

    fprintf('  %s — top-3 by mean rank: %s, %s, %s\n', metric_labels{m}, ...
        T_C.Family{1}, T_C.Family{2}, T_C.Family{3});
end

% ════════════════════════════════════════════════════════════════════════════
%  TABLE D — Edge persistence
% ════════════════════════════════════════════════════════════════════════════
fprintf('[robustness] Building Table D (edge persistence)...\n');
[ii_all, jj_all] = find(triu(edge_count > 0));
n_ep = numel(ii_all);

famA_ep  = families(ii_all);
famB_ep  = families(jj_all);
cnt_ep   = edge_count(sub2ind([N,N], ii_all, jj_all));
mean_ep  = edge_sumDV(sub2ind([N,N], ii_all, jj_all)) ./ cnt_ep;
% std = sqrt(E[X²] - E[X]²), guarded for n=1
var_ep   = edge_sumDV2(sub2ind([N,N], ii_all, jj_all)) ./ cnt_ep - mean_ep.^2;
std_ep   = sqrt(max(var_ep, 0));
frac_ep  = cnt_ep / nValid;

T_D = table(famA_ep(:), famB_ep(:), cnt_ep(:), frac_ep(:), mean_ep(:), std_ep(:), ...
    'VariableNames', {'FamilyA','FamilyB','N_runs','FracRuns','MeanDV_mps','StdDV_mps'});
T_D = sortrows(T_D, 'N_runs', 'descend');
writetable(T_D, fullfile(outDir, 'TableD_edge_persistence.csv'));
fprintf('[robustness] %d unique edges found across all runs.\n\n', n_ep);

% ════════════════════════════════════════════════════════════════════════════
%  FIGURE 1 — Rank heatmap
% ════════════════════════════════════════════════════════════════════════════
fprintf('[robustness] Generating figures...\n');

fig1 = figure('Name','Rank Heatmap','NumberTitle','off', ...
    'Units','normalized','Position',[0.02 0.05 0.95 0.85]);
run_short = cellfun(@(s) strrep(strrep(s,'Run','R'),'_',' '), run_ids(valid_idx), ...
    'UniformOutput', false);

for m = 1:nMetrics
    subplot(2, 3, m);
    rnk_v = rank_all{m}(:, valid_idx);
    imagesc(rnk_v);
    colormap(gca, flipud(parula(N)));
    clim([1 N]);
    cb = colorbar;
    cb.Label.String = 'Rank (1=best)';
    set(gca, 'XTick', 1:nValid, 'XTickLabel', run_short, 'XTickLabelRotation', 45, ...
             'YTick', 1:N, 'YTickLabel', short_names, 'FontSize', 7);
    title(metric_labels{m}, 'FontSize', 9, 'FontWeight', 'bold');
    xlabel('Run'); ylabel('Family');
end

sgtitle('Family Rank per Metric Across Verification Runs', 'FontWeight','bold');

if SAVE_FIGS
    saveas(fig1, fullfile(outDir, 'Fig1_rank_heatmap.png'));
    saveas(fig1, fullfile(outDir, 'Fig1_rank_heatmap.fig'));
end

% ════════════════════════════════════════════════════════════════════════════
%  FIGURE 2 — Top-5 bar charts with std-rank error bars
% ════════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Top-5 Robustness','NumberTitle','off', ...
    'Units','normalized','Position',[0.02 0.05 0.95 0.85]);

for m = 1:nMetrics
    rnk_v = rank_all{m}(:, valid_idx);
    m_r   = nanmean(rnk_v, 2);
    s_r   = nanstd(rnk_v, 0, 2);
    [m_r_s, ord] = sort(m_r, 'ascend');
    s_r_s = s_r(ord);
    top5  = min(5, N);
    idx5  = ord(1:top5);

    subplot(2, 3, m);
    bar_h = bar(1:top5, m_r_s(1:top5), 0.6, 'FaceColor', [0.3 0.6 0.9]);
    hold on;
    errorbar(1:top5, m_r_s(1:top5), s_r_s(1:top5), 'k.', 'LineWidth', 1.2);
    hold off;
    set(gca, 'XTick', 1:top5, 'XTickLabel', short_names(idx5), ...
             'XTickLabelRotation', 30, 'FontSize', 7);
    ylabel('Mean Rank (↓ better)');
    title(metric_labels{m}, 'FontSize', 9, 'FontWeight', 'bold');
    ylim([0 N]);
    grid on;
end

sgtitle('Top-5 Families by Mean Rank (error bars = std across runs)', 'FontWeight','bold');

if SAVE_FIGS
    saveas(fig2, fullfile(outDir, 'Fig2_top5_robustness.png'));
    saveas(fig2, fullfile(outDir, 'Fig2_top5_robustness.fig'));
end

% ════════════════════════════════════════════════════════════════════════════
%  FIGURE 3 — Edge-persistence heatmap
% ════════════════════════════════════════════════════════════════════════════
persist_frac = edge_count / nValid;   % fraction of runs [0,1]; symmetric
persist_sym  = persist_frac + persist_frac';

fig3 = figure('Name','Edge Persistence','NumberTitle','off', ...
    'Units','normalized','Position',[0.1 0.1 0.55 0.55]);
imagesc(persist_sym);
colormap(flipud(hot));
clim([0 2]);   % symmetric so max is 2 (both [i,j] and [j,i] counted)
cb = colorbar;
cb.Label.String = 'Fraction of runs (×2 for symmetry)';
set(gca, 'XTick', 1:N, 'XTickLabel', short_names, 'XTickLabelRotation', 45, ...
         'YTick', 1:N, 'YTickLabel', short_names, 'FontSize', 8);
title(sprintf('Edge Persistence Across %d Runs  (bright = always present)', nValid), ...
    'FontWeight', 'bold');
axis square;

if SAVE_FIGS
    saveas(fig3, fullfile(outDir, 'Fig3_edge_persistence.png'));
    saveas(fig3, fullfile(outDir, 'Fig3_edge_persistence.fig'));
end

% ════════════════════════════════════════════════════════════════════════════
%  FIGURE 4 — Baseline network (Run 0 if available, else first valid)
% ════════════════════════════════════════════════════════════════════════════
base_r = valid_idx(1);   % use first valid run as baseline
DV_base = DV_all{base_r};
if ~isempty(DV_base)
    A_base  = isfinite(DV_base) & (DV_base <= DV_budget) & ~eye(N);
    [ii_b, jj_b] = find(triu(A_base, 1));
    dv_b = DV_base(sub2ind([N,N], ii_b, jj_b));

    fig4 = figure('Name', sprintf('Baseline Network (%s)', run_ids{base_r}), ...
        'NumberTitle','off','Units','normalized','Position',[0.1 0.1 0.55 0.55]);

    % Circular layout
    theta_pos = linspace(0, 2*pi, N+1)';
    theta_pos = theta_pos(1:N);
    xn = cos(theta_pos);
    yn = sin(theta_pos);

    hold on;
    dv_max = max(dv_b);
    dv_min = min(dv_b);
    dv_range = max(dv_max - dv_min, 1);
    for e = 1:numel(ii_b)
        lw = 4 * (1 - (dv_b(e) - dv_min) / dv_range) + 0.5;
        col = [0.3 + 0.6*(dv_b(e)-dv_min)/dv_range, ...
               0.7 - 0.5*(dv_b(e)-dv_min)/dv_range, ...
               0.9 - 0.6*(dv_b(e)-dv_min)/dv_range];
        plot([xn(ii_b(e)), xn(jj_b(e))], [yn(ii_b(e)), yn(jj_b(e))], ...
            '-', 'LineWidth', lw, 'Color', col);
    end
    scatter(xn, yn, 120, 'k', 'filled');
    for i = 1:N
        text(xn(i)*1.13, yn(i)*1.13, short_names{i}, ...
            'HorizontalAlignment','center', 'FontSize', 7, 'FontWeight','bold');
    end
    hold off;
    axis equal off;
    title(sprintf('Baseline Network — %s\n%d edges, DV \\leq %.0f m/s', ...
        run_ids{base_r}, numel(ii_b), DV_budget), 'FontWeight','bold');
    % Colourbar legend: thin=high DV, thick=low DV
    annotation('textbox',[0.75 0.07 0.2 0.08],'String', ...
        'Thick/blue = low DV\nThin/red = high DV', ...
        'FitBoxToText','on','EdgeColor','none','FontSize',7);

    if SAVE_FIGS
        saveas(fig4, fullfile(outDir, 'Fig4_baseline_network.png'));
        saveas(fig4, fullfile(outDir, 'Fig4_baseline_network.fig'));
    end
end

% ════════════════════════════════════════════════════════════════════════════
%  FIGURE 5 — Coarse vs fine comparison
%  Panels: Run 0 (baseline) | Run 5 (combined coarse) | Run 6 (combined fine)
% ════════════════════════════════════════════════════════════════════════════
compare_ids = {'Run0_Baseline', 'Run5_CombinedCoarse', 'Run6_CombinedFine'};
compare_titles = {'Run 0  Baseline', 'Run 5  Combined Coarse', 'Run 6  Combined Fine'};

compare_ridx = cellfun(@(id) find(strcmp(run_ids, id), 1), compare_ids, ...
    'UniformOutput', false);
has_compare = ~cellfun(@isempty, compare_ridx);

if any(has_compare)
    fig5 = figure('Name','Coarse vs Fine Comparison','NumberTitle','off', ...
        'Units','normalized','Position',[0.02 0.1 0.95 0.45]);

    theta_pos = linspace(0, 2*pi, N+1)';
    theta_pos = theta_pos(1:N);
    xn = cos(theta_pos);
    yn = sin(theta_pos);

    for panel = 1:3
        if ~has_compare(panel), continue; end
        ri = compare_ridx{panel};
        DV_p = DV_all{ri};
        if isempty(DV_p), continue; end

        A_p = isfinite(DV_p) & (DV_p <= DV_budget) & ~eye(N);
        [ii_p, jj_p] = find(triu(A_p, 1));
        dv_p = DV_p(sub2ind([N,N], ii_p, jj_p));

        subplot(1, 3, panel);
        hold on;
        if ~isempty(dv_p)
            dv_lo = min(dv_p);  dv_hi = max(dv_p);  dv_rng = max(dv_hi-dv_lo, 1);
            for e = 1:numel(ii_p)
                lw = 3*(1-(dv_p(e)-dv_lo)/dv_rng) + 0.5;
                col = [0.3+0.6*(dv_p(e)-dv_lo)/dv_rng, ...
                       0.7-0.5*(dv_p(e)-dv_lo)/dv_rng, ...
                       0.9-0.6*(dv_p(e)-dv_lo)/dv_rng];
                plot([xn(ii_p(e)), xn(jj_p(e))], [yn(ii_p(e)), yn(jj_p(e))], ...
                    '-', 'LineWidth', lw, 'Color', col);
            end
        end
        scatter(xn, yn, 80, 'k', 'filled');
        for i = 1:N
            text(xn(i)*1.18, yn(i)*1.18, short_names{i}, ...
                'HorizontalAlignment','center','FontSize',6,'FontWeight','bold');
        end
        hold off;
        axis equal off;
        title(sprintf('%s\n(%d edges)', compare_titles{panel}, numel(ii_p)), ...
            'FontSize', 8, 'FontWeight', 'bold');
    end

    sgtitle(sprintf('Coarse vs Fine Network Comparison  (DV \\leq %.0f m/s)', DV_budget), ...
        'FontWeight', 'bold');

    if SAVE_FIGS
        saveas(fig5, fullfile(outDir, 'Fig5_coarse_fine_comparison.png'));
        saveas(fig5, fullfile(outDir, 'Fig5_coarse_fine_comparison.fig'));
    end
end

% ════════════════════════════════════════════════════════════════════════════
%  DONE
% ════════════════════════════════════════════════════════════════════════════
fprintf('\n[robustness] ══════════ ANALYSIS COMPLETE ══════════\n');
fprintf('[robustness] Output dir: %s\n', outDir);
fprintf('[robustness] Files written:\n');
fprintf('   TableA_graph_summary.csv\n');
for m = 1:nMetrics
    fprintf('   TableB_%s_per_run.csv\n', metric_ids{m});
    fprintf('   TableC_%s_robustness.csv\n', metric_ids{m});
end
fprintf('   TableD_edge_persistence.csv\n');
if SAVE_FIGS
    fprintf('   Fig1_rank_heatmap.{png,fig}\n');
    fprintf('   Fig2_top5_robustness.{png,fig}\n');
    fprintf('   Fig3_edge_persistence.{png,fig}\n');
    fprintf('   Fig4_baseline_network.{png,fig}\n');
    fprintf('   Fig5_coarse_fine_comparison.{png,fig}\n');
end

% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FALLBACK FUNCTIONS  (used only if net_* functions are not on path)
% ════════════════════════════════════════════════════════════════════════════

function dist = local_floyd_warshall(W, N)
dist = W;
for k = 1:N
    dist = min(dist, dist(:,k) + dist(k,:));
end
dist(1:N+1:end) = 0;
end

function lcc_sz = local_lcc(A, N)
visited = false(N, 1);
best = 0;
for start = 1:N
    if visited(start), continue; end
    queue = start;
    comp  = 0;
    while ~isempty(queue)
        v = queue(1);  queue = queue(2:end);
        if visited(v), continue; end
        visited(v) = true;
        comp = comp + 1;
        nbrs = find(A(v,:));
        queue = [queue, nbrs(~visited(nbrs))]; %#ok<AGROW>
    end
    best = max(best, comp);
end
lcc_sz = best;
end

function btw = local_betweenness(A, dist, N)
% Simplified node betweenness: fraction of (s,t) shortest paths through v.
btw = zeros(N, 1);
n_paths = 0;
for s = 1:N
    for t = s+1:N
        if ~isfinite(dist(s,t)), continue; end
        n_paths = n_paths + 1;
        for v = 1:N
            if v == s || v == t, continue; end
            if isfinite(dist(s,v)) && isfinite(dist(v,t)) && ...
               abs(dist(s,v) + dist(v,t) - dist(s,t)) < 1e-9 * dist(s,t)
                btw(v) = btw(v) + 1;
            end
        end
    end
end
if n_paths > 0
    btw = btw / n_paths;
end
end
