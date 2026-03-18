function net_plot_winner_map(plot_type, map_data, DV_vec, Tmax_vec, ...
    lcc_full_map, budget_pairs_map, short_names, N_fam, out_dir)
%NET_PLOT_WINNER_MAP  Create, save (PDF + PNG + .fig) one network-centrality figure.
%
%   Saves <out_dir>/<plot_type>.{pdf,png,fig}.
%
% Inputs
%   plot_type        char   one of:
%                      'edges_count_map'
%                      'budget_feasible_pairs_map'
%                      'lcc_map'
%                      'articulation_map'
%                      'winner_map_strength'
%                      'winner_map_harmonic_closeness'
%                      'winner_map_betweenness'
%                      'strength_contour'
%                      'harmonic_closeness_contour'
%   map_data         struct  fields depend on plot_type (see switch below)
%   DV_vec           [nDV×1]   x-axis (DVcap_true_mps, m/s)
%   Tmax_vec         [nTmax×1] y-axis (Tmax_true_days, days)
%   lcc_full_map     [nDV×nTmax] binary; 1 where entire graph is connected
%   budget_pairs_map [nDV×nTmax] undirected budget-feasible pair count (max 78)
%   short_names      {N_fam×1}  short family labels
%   N_fam            scalar     number of families (13)
%   out_dir          char       output directory (created if absent)

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

nDV   = numel(DV_vec);
nTmax = numel(Tmax_vec);

% Cell half-widths — used for tie-cell patch rendering
hw = i_cell_hw(DV_vec);
hh = i_cell_hw(Tmax_vec);

fam_colors = i_fam_colors(N_fam);

fig = figure('Visible', 'off', 'Color', 'white', ...
             'Units', 'pixels', 'Position', [100 100 860 620]);
ax  = axes(fig);

% ── MAX PAIRS THRESHOLD (upper triangle of 13×13 network) ────────────────────
MAX_PAIRS = N_fam * (N_fam - 1) / 2;   % = 78 for N_fam=13

switch plot_type

    % ═══════════════════════════════════════════════════════════════════════
    %  COUNT HEATMAPS (direct edges  /  budget-feasible paths)
    % ═══════════════════════════════════════════════════════════════════════
    case {'edges_count_map', 'budget_feasible_pairs_map'}

        if strcmp(plot_type, 'edges_count_map')
            count_data = map_data.direct_pairs;   % [nDV×nTmax] undirected, max 78
            cb_label   = 'Direct family-pair count';
            ttl        = 'Direct family-pair count';
        else
            count_data = map_data.budget_pairs;   % [nDV×nTmax] max 78
            cb_label   = 'Budget-feasible transfer pairs';
            ttl        = 'Budget-feasible transfer pairs';
        end

        skip_mask = map_data.skip_map;
        Z = double(count_data);
        Z(skip_mask) = NaN;
        Z = Z';   % [nTmax × nDV] for imagesc

        imagesc(ax, DV_vec, Tmax_vec, Z);
        axis(ax, 'xy');
        set(ax, 'Color', [1 1 1]);   % white for NaN/skip cells
        colormap(ax, parula(256));
        cb = colorbar(ax);
        cb.Label.String = cb_label;
        cb.Label.FontSize = 9;

        hold(ax, 'on');

        % Solid black labeled contours at 10, 20, … up to MAX_PAIRS-8
        max_level = floor((MAX_PAIRS - 1) / 10) * 10;   % = 70 for MAX_PAIRS=78
        lvls = 10 : 10 : max_level;
        if ~isempty(lvls) && any(isfinite(Z(:)))
            [~, hc_cnt] = contour(ax, DV_vec, Tmax_vec, Z, lvls, 'k-', 'LineWidth', 0.8);
            hc_cnt.ShowText    = 'on';
            hc_cnt.LabelSpacing = 500;
        end

        % Two special contours (black dashed LCC=13, red dashed budget=78)
        [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                lcc_full_map, budget_pairs_map, MAX_PAIRS);

        i_add_contour_legend(ax, h_lcc, h_bp);
        title(ax, ttl, 'FontSize', 11);

    % ═══════════════════════════════════════════════════════════════════════
    %  WINNER MAPS
    % ═══════════════════════════════════════════════════════════════════════
    case {'winner_map_strength', 'winner_map_harmonic_closeness', ...
          'winner_map_betweenness'}

        winner_idx   = map_data.winner_idx;    % [nDV×nTmax]  1-N, 0=tie, NaN=skip
        winner_names = map_data.winner_names;  % {nDV×nTmax}  semicolon-joined if tie

        % Build RGB background: winners → family colour; ties/skip → white
        rgb = ones(nTmax, nDV, 3);
        for di = 1:nDV
            for dj = 1:nTmax
                idx = winner_idx(di, dj);
                if isfinite(idx) && idx > 0
                    rgb(dj, di, :) = fam_colors(idx, :);
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');
        hold(ax, 'on');

        % Overlay tie cells as diagonal-split (2-way) or pie-slice (3+-way) patches
        i_draw_tie_cells(ax, winner_idx, winner_names, short_names, fam_colors, ...
                         DV_vec, Tmax_vec, hw, hh);

        % Two special contours
        [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                lcc_full_map, budget_pairs_map, MAX_PAIRS);

        % Filtered legend: only families that actually win somewhere
        present_idx = unique(winner_idx(isfinite(winner_idx) & winner_idx > 0));
        leg_items   = {};
        for k = 1:numel(present_idx)
            p = patch(ax, NaN, NaN, fam_colors(present_idx(k), :), ...
                      'EdgeColor', 'none', ...
                      'DisplayName', short_names{present_idx(k)});
            leg_items{end+1} = p; %#ok<AGROW>
        end
        if ~isempty(h_lcc), leg_items{end+1} = h_lcc; end
        if ~isempty(h_bp),  leg_items{end+1} = h_bp;  end
        if ~isempty(leg_items)
            legend(ax, [leg_items{:}], 'Location', 'eastoutside', ...
                   'FontSize', 8, 'Box', 'on');
        end

        switch plot_type
            case 'winner_map_strength'
                ttl = 'Strength winner';
            case 'winner_map_harmonic_closeness'
                ttl = 'Harmonic closeness winner';
            case 'winner_map_betweenness'
                ttl = 'Betweenness winner';
        end
        title(ax, ttl, 'FontSize', 11);

    % ═══════════════════════════════════════════════════════════════════════
    %  LCC SIZE MAP
    % ═══════════════════════════════════════════════════════════════════════
    case 'lcc_map'

        lcc_size  = map_data.lcc_size;   % [nDV×nTmax]
        skip_mask = (lcc_size == 0);

        valid_sz = lcc_size(~skip_mask);
        if isempty(valid_sz)
            min_sz = 1;  max_sz = 1;
        else
            min_sz = min(valid_sz(:));
            max_sz = max(valid_sz(:));
        end

        % Float matrix (NaN for skip)
        Z = double(lcc_size);
        Z(skip_mask) = NaN;
        Z = Z';   % [nTmax × nDV]

        % Discrete colormap (parula slice)
        full_cmap = i_lcc_colors(N_fam);
        sub_cmap  = full_cmap(min_sz:max_sz, :);

        imagesc(ax, DV_vec, Tmax_vec, Z, [min_sz - 0.5, max_sz + 0.5]);
        axis(ax, 'xy');
        set(ax, 'Color', [1 1 1]);
        colormap(ax, sub_cmap);

        cb = colorbar(ax);
        cb.Ticks      = min_sz : max_sz;
        cb.TickLabels = arrayfun(@num2str, min_sz:max_sz, 'UniformOutput', false);
        cb.Label.String   = 'Largest connected component size';
        cb.Label.FontSize = 9;

        hold(ax, 'on');
        [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                lcc_full_map, budget_pairs_map, MAX_PAIRS);
        i_add_contour_legend(ax, h_lcc, h_bp);
        title(ax, 'Largest connected component size', 'FontSize', 11);

    % ═══════════════════════════════════════════════════════════════════════
    %  ARTICULATION-POINT COUNT MAP
    % ═══════════════════════════════════════════════════════════════════════
    case 'articulation_map'

        ap_count = double(map_data.ap_count);
        skip_mask = map_data.skip_map;

        max_ap   = max(ap_count(:));
        if max_ap == 0, max_ap = 1; end

        n_colors  = max_ap + 1;
        ap_cmap   = parula(max(n_colors, 2));

        rgb = ones(nTmax, nDV, 3);
        for di = 1:nDV
            for dj = 1:nTmax
                if skip_mask(di, dj)
                    rgb(dj, di, :) = [1 1 1];
                else
                    ci = ap_count(di, dj) + 1;
                    ci = max(1, min(ci, n_colors));
                    rgb(dj, di, :) = ap_cmap(ci, :);
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');
        hold(ax, 'on');
        [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                lcc_full_map, budget_pairs_map, MAX_PAIRS);

        n_show = min(n_colors, 8);
        hp = gobjects(n_show, 1);
        for k = 1:n_show
            ci = round((k-1) * (n_colors-1) / max(n_show-1,1)) + 1;
            ci = max(1, min(ci, n_colors));
            hp(k) = patch(ax, NaN, NaN, ap_cmap(ci,:), 'EdgeColor', 'none', ...
                          'DisplayName', sprintf('%d', ci-1));
        end
        leg_items = num2cell(hp);
        if ~isempty(h_lcc), leg_items{end+1} = h_lcc; end
        if ~isempty(h_bp),  leg_items{end+1} = h_bp;  end
        lgd = legend(ax, [leg_items{:}], 'Location', 'eastoutside', ...
                     'FontSize', 7, 'Box', 'on');
        lgd.Title.String = '# articulation pts';
        title(ax, 'Articulation-point count per snapshot', 'FontSize', 11);

    % ═══════════════════════════════════════════════════════════════════════
    %  CONTINUOUS METRIC CONTOUR MAPS
    % ═══════════════════════════════════════════════════════════════════════
    case {'strength_contour', 'harmonic_closeness_contour'}

        vals_3d  = map_data.values;    % [N × nDV × nTmax]
        skip_mask = map_data.skip_map;

        max_vals = squeeze(max(vals_3d, [], 1));   % [nDV × nTmax]
        max_vals(skip_mask) = NaN;
        Z = max_vals';   % [nTmax × nDV]

        imagesc(ax, DV_vec, Tmax_vec, Z);
        axis(ax, 'xy');
        set(ax, 'Color', [1 1 1]);
        colormap(ax, parula(256));
        cb = colorbar(ax);

        switch plot_type
            case 'strength_contour'
                ttl    = 'Strength (max over families)';
                cb_lbl = 'Strength';
            case 'harmonic_closeness_contour'
                ttl    = 'Harmonic closeness (max over families)';
                cb_lbl = 'Harmonic closeness';
        end
        cb.Label.String   = cb_lbl;
        cb.Label.FontSize = 9;

        hold(ax, 'on');
        valid_data = Z(isfinite(Z));
        if ~isempty(valid_data)
            vmin = min(valid_data);  vmax = max(valid_data);
            if vmax > vmin
                step = i_nice_step((vmax - vmin) / 8);
                lvls = (ceil(vmin/step)*step) : step : vmax;
                if numel(lvls) > 1
                    [~, hc] = contour(ax, DV_vec, Tmax_vec, Z, lvls, 'k-', 'LineWidth', 0.8);
                    hc.ShowText    = 'on';
                    hc.LabelSpacing = 500;
                end
            end
        end

        [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                lcc_full_map, budget_pairs_map, MAX_PAIRS);
        i_add_contour_legend(ax, h_lcc, h_bp);
        title(ax, ttl, 'FontSize', 11);

    otherwise
        error('net_plot_winner_map: unknown plot_type ''%s''.', plot_type);
end

% ── Common axis labels + formatting ─────────────────────────────────────────
xlabel(ax, '\DeltaV_{cap}  [m/s]',  'FontSize', 10);
ylabel(ax, 'T_{cap}  [days]',       'FontSize', 10);
set(ax, 'FontSize', 9, 'Box', 'on', 'TickDir', 'out', 'Layer', 'top');
xlim(ax, [DV_vec(1), DV_vec(end)]);
ylim(ax, [Tmax_vec(1), Tmax_vec(end)]);

% ── Save: PDF + PNG + .fig ───────────────────────────────────────────────────
base = fullfile(out_dir, plot_type);

set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8.6 6.2], ...
         'PaperSize', [8.6 6.2]);
print(fig, [base '.pdf'], '-dpdf', '-bestfit');
print(fig, [base '.png'], '-dpng', '-r150');
saveas(fig, [base '.fig']);

close(fig);

end   % end main function

% ═══════════════════════════════════════════════════════════════════════════
%  Local helpers
% ═══════════════════════════════════════════════════════════════════════════

function [h_lcc, h_bp] = i_overlay_two_contours(ax, DV_vec, Tmax_vec, ...
                                                  lcc_full_map, budget_pairs_map, max_pairs)
%I_OVERLAY_TWO_CONTOURS  Draw two special contour lines on ax:
%   1. Black dashed, inline-labelled "LCC=13" — where the network first becomes
%      fully connected (lcc_full_map transitions 0→1).
%   2. Red dashed, inline-labelled with the max-pairs value — where all
%      budget-feasible unordered pairs first reach the maximum (max_pairs, e.g. 78).

hold(ax, 'on');
h_lcc = [];
h_bp  = [];

% ── 1. Black dashed: LCC fully connected ──────────────────────────────────────
Z_lcc = lcc_full_map';   % [nTmax × nDV]
if any(Z_lcc(:) == 0) && any(Z_lcc(:) == 1)
    [C1, h1] = contour(ax, DV_vec, Tmax_vec, Z_lcc, [0.5 0.5], ...
                       'k--', 'LineWidth', 1.5);
    cl1 = clabel(C1, h1, 'FontSize', 7, 'Color', 'k');
    for ii = 1:numel(cl1)
        if isgraphics(cl1(ii), 'text')
            set(cl1(ii), 'String', 'LCC=13', 'BackgroundColor', 'none');
        end
    end
    % Dummy line for legend
    h_lcc = plot(ax, NaN, NaN, 'k--', 'LineWidth', 1.5, 'DisplayName', 'LCC=13');
end

% ── 2. Red dashed: all pairs budget-reachable ─────────────────────────────────
if ~isempty(budget_pairs_map) && any(isfinite(budget_pairs_map(:)))
    Z_bp = budget_pairs_map';   % [nTmax × nDV]
    thresh = max_pairs - 0.5;
    if any(Z_bp(:) < thresh) && any(Z_bp(:) >= thresh)
        [C2, h2] = contour(ax, DV_vec, Tmax_vec, Z_bp, [thresh thresh], ...
                           'r--', 'LineWidth', 2.0);
        lbl_str = sprintf('%d', max_pairs);
        cl2 = clabel(C2, h2, 'FontSize', 7, 'Color', 'r');
        for ii = 1:numel(cl2)
            if isgraphics(cl2(ii), 'text')
                set(cl2(ii), 'String', lbl_str, 'BackgroundColor', 'none');
            end
        end
        h_bp = plot(ax, NaN, NaN, 'r--', 'LineWidth', 2.0, ...
                    'DisplayName', sprintf('Budget: %d pairs', max_pairs));
    end
end

end   % end i_overlay_two_contours

% ─────────────────────────────────────────────────────────────────────────────

function i_add_contour_legend(ax, h_lcc, h_bp)
%I_ADD_CONTOUR_LEGEND  Add a compact legend showing only the two special contours.
items = {};
if ~isempty(h_lcc), items{end+1} = h_lcc; end
if ~isempty(h_bp),  items{end+1} = h_bp;  end
if ~isempty(items)
    legend(ax, [items{:}], 'Location', 'best', 'FontSize', 7, 'Box', 'on');
end
end

% ─────────────────────────────────────────────────────────────────────────────

function i_draw_tie_cells(ax, winner_idx, winner_names, short_names, fam_colors, ...
                           DV_vec, Tmax_vec, hw, hh)
%I_DRAW_TIE_CELLS  Overlay diagonal-split (2-way) or pie (3+-way) patches on ax.

for di = 1:numel(DV_vec)
    for dj = 1:numel(Tmax_vec)
        if winner_idx(di, dj) ~= 0, continue; end   % only tie cells (0); NaN~=0 is fine

        name_str = winner_names{di, dj};
        if isempty(name_str), continue; end

        parts   = strsplit(name_str, ';');
        n_tied  = numel(parts);
        tidx    = zeros(1, n_tied);
        for f = 1:n_tied
            idx = find(strcmp(short_names, strtrim(parts{f})), 1);
            tidx(f) = idx;
        end
        valid = tidx > 0 & tidx <= size(fam_colors, 1);
        if ~any(valid), continue; end
        tidx   = tidx(valid);
        n_tied = numel(tidx);

        cx = DV_vec(di);
        cy = Tmax_vec(dj);

        if n_tied == 2
            % Split along BL→TR diagonal
            % Upper-left triangle: BL, TL, TR
            patch(ax, [cx-hw  cx-hw  cx+hw], [cy-hh  cy+hh  cy+hh], ...
                  fam_colors(tidx(1),:), 'EdgeColor', 'none');
            % Lower-right triangle: BL, BR, TR
            patch(ax, [cx-hw  cx+hw  cx+hw], [cy-hh  cy-hh  cy+hh], ...
                  fam_colors(tidx(2),:), 'EdgeColor', 'none');

        else
            % Equal pie slices
            n_pts = 24;
            r     = min(hw, hh) * 0.95;
            for s = 1:n_tied
                th1 = (s-1) * 2*pi/n_tied - pi/2;
                th2 =  s    * 2*pi/n_tied - pi/2;
                th  = linspace(th1, th2, n_pts);
                xv  = [cx,  cx + r*cos(th)];
                yv  = [cy,  cy + r*sin(th)];
                patch(ax, xv, yv, fam_colors(tidx(s),:), 'EdgeColor', 'none');
            end
        end
    end
end
end   % end i_draw_tie_cells

% ─────────────────────────────────────────────────────────────────────────────

function hw = i_cell_hw(vec)
%I_CELL_HW  Half-width of one grid cell along vec.
if numel(vec) > 1
    hw = abs(vec(2) - vec(1)) / 2;
else
    hw = max(abs(vec(1)) * 0.05, 0.5);
end
end

% ─────────────────────────────────────────────────────────────────────────────

function step = i_nice_step(raw_step)
%I_NICE_STEP  Round raw_step up to a "nice" number (1, 2, 5 × 10^k).
if raw_step <= 0, step = 1; return; end
mag  = 10^floor(log10(raw_step));
frac = raw_step / mag;
if     frac <= 1, step = 1 * mag;
elseif frac <= 2, step = 2 * mag;
elseif frac <= 5, step = 5 * mag;
else,             step = 10 * mag;
end
end

% ─────────────────────────────────────────────────────────────────────────────

function c = i_fam_colors(N)
%I_FAM_COLORS  N visually distinct colours for the N families.
base = lines(max(N, 7));
if N <= size(base, 1)
    c = base(1:N, :);
else
    c = hsv(N);
end
end

% ─────────────────────────────────────────────────────────────────────────────

function c = i_lcc_colors(N)
%I_LCC_COLORS  N sequential colours (parula) for LCC-size values 1..N.
c = parula(max(N, 2));
end
