function net_plot_winner_map(plot_type, map_data, DV_vec, Tmax_vec, ...
    lcc_full_map, short_names, N_fam, out_dir)
%NET_PLOT_WINNER_MAP  Create, save (PDF + PNG + .fig) one network-centrality figure.
%
%   Called once per figure type.  Saves <out_dir>/<plot_type>.{pdf,png,fig}.
%
% Inputs
%   plot_type     char    one of:
%                   'edges_count_map'
%                   'lcc_map'
%                   'articulation_map'
%                   'winner_map_strength'
%                   'winner_map_harmonic_closeness'
%                   'winner_map_betweenness'
%                   'strength_contour'
%                   'harmonic_closeness_contour'
%                   'betweenness_contour'
%   map_data      struct  fields depend on plot_type (see below)
%   DV_vec        [nDV×1]   x-axis values (DVcap_true_mps, m/s)
%   Tmax_vec      [nTmax×1] y-axis values (Tmax_true_days, days)
%   lcc_full_map  [nDV×nTmax] binary; 1 where entire graph is connected
%   short_names   {N_fam×1} short family labels
%   N_fam         scalar    number of families (13)
%   out_dir       char      output directory (created if absent)
%
% map_data fields by plot_type
%   edges_count_map          .edges_kept [nDV×nTmax]  .skip_map [nDV×nTmax]
%   lcc_map                  .lcc_size   [nDV×nTmax]
%   articulation_map         .ap_count   [nDV×nTmax]  .skip_map [nDV×nTmax]
%   winner_map_*             .winner_idx [nDV×nTmax] (1-N, 0=Tie, NaN=skip)
%                            .tie_size   [nDV×nTmax]
%   *_contour                .values     [N×nDV×nTmax]  .skip_map [nDV×nTmax]

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

nDV   = numel(DV_vec);
nTmax = numel(Tmax_vec);

% ── Build figure ─────────────────────────────────────────────────────────────
fig = figure('Visible', 'off', 'Color', 'white', ...
             'Units', 'pixels', 'Position', [100 100 820 600]);
ax  = axes(fig);

switch plot_type

    % ── Direct family-pair count ──────────────────────────────────────────
    case 'edges_count_map'

        edge_data = double(map_data.edges_kept);   % [nDV×nTmax]
        skip_map  = map_data.skip_map;

        % Replace skip cells with NaN for display
        edge_data(skip_map) = NaN;

        Z = edge_data';   % [nTmax × nDV] for imagesc

        imagesc(ax, DV_vec, Tmax_vec, Z);
        axis(ax, 'xy');
        colormap(ax, parula(256));
        cb = colorbar(ax);
        cb.Label.String = 'Direct family-pair count';

        % Contour overlay
        hold(ax, 'on');
        valid_min = min(edge_data(~skip_map));
        valid_max = max(edge_data(~skip_map));
        if ~isempty(valid_min) && valid_max > valid_min
            step = max(1, round((valid_max - valid_min) / 8));
            levels = (ceil(valid_min/step)*step) : step : valid_max;
            if numel(levels) > 1
                [~, hc] = contour(ax, DV_vec, Tmax_vec, Z, levels, 'k-');
                hc.ShowText = 'on';
                hc.LabelSpacing = 400;
            end
        end

        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);
        title(ax, 'Direct family-pair count', 'FontSize', 11);

    % ── Winner maps (strength / harmonic_closeness / betweenness) ────────
    case {'winner_map_strength', 'winner_map_harmonic_closeness', ...
          'winner_map_betweenness'}

        winner_idx = map_data.winner_idx;   % [nDV×nTmax]

        % 13 distinct family colours (lines(13) palette)
        fam_colors = i_fam_colors(N_fam);
        tie_color  = [0.65 0.65 0.65];      % grey

        % Build RGB image [nTmax × nDV × 3]
        rgb = ones(nTmax, nDV, 3);           % default white = skip

        for di = 1:nDV
            for dj = 1:nTmax
                idx = winner_idx(di, dj);
                if isnan(idx)
                    rgb(dj, di, :) = [1 1 1];
                elseif idx == 0
                    rgb(dj, di, :) = tie_color;
                else
                    rgb(dj, di, :) = fam_colors(idx, :);
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');

        % Overlay LCC full contour (red dashed)
        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

        % For betweenness: overlay hollow circles at tie cells
        if strcmp(plot_type, 'winner_map_betweenness')
            hold(ax, 'on');
            [ti, tj] = find(winner_idx == 0);
            if ~isempty(ti)
                scatter(ax, DV_vec(ti), Tmax_vec(tj), 60, 'wo', ...
                        'LineWidth', 1.5, 'MarkerFaceColor', 'none');
            end
        end

        % Legend: coloured patches
        hold(ax, 'on');
        hp = gobjects(N_fam + 2, 1);
        for k = 1:N_fam
            hp(k) = patch(ax, NaN, NaN, fam_colors(k,:), ...
                          'EdgeColor', 'none', 'DisplayName', short_names{k});
        end
        hp(N_fam+1) = patch(ax, NaN, NaN, tie_color, ...
                            'EdgeColor', 'none', 'DisplayName', 'Tie');
        hp(N_fam+2) = patch(ax, NaN, NaN, [1 1 1], ...
                            'EdgeColor', [0.8 0.8 0.8], 'DisplayName', 'Skip');

        legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'on');

        switch plot_type
            case 'winner_map_strength'
                ttl = 'Strength winner';
            case 'winner_map_harmonic_closeness'
                ttl = 'Harmonic closeness winner';
            case 'winner_map_betweenness'
                ttl = 'Betweenness winner';
        end
        title(ax, ttl, 'FontSize', 11);

    % ── Continuous metric heatmaps (strength / HC / betweenness) ─────────
    case {'strength_contour', 'harmonic_closeness_contour', 'betweenness_contour'}

        vals_3d  = map_data.values;    % [N × nDV × nTmax]
        skip_map = map_data.skip_map;  % [nDV × nTmax]

        % Max-over-families for the heatmap (shows the "best" value at each budget)
        max_vals = squeeze(max(vals_3d, [], 1));   % [nDV × nTmax]
        max_vals(skip_map) = NaN;

        Z = max_vals';   % [nTmax × nDV]

        imagesc(ax, DV_vec, Tmax_vec, Z);
        axis(ax, 'xy');
        colormap(ax, parula(256));
        cb = colorbar(ax);

        switch plot_type
            case 'strength_contour'
                ttl     = 'Strength (max over families)';
                cb_lbl  = 'Strength';
            case 'harmonic_closeness_contour'
                ttl     = 'Harmonic closeness (max over families)';
                cb_lbl  = 'Harmonic closeness';
            case 'betweenness_contour'
                ttl     = 'Betweenness (max over families)';
                cb_lbl  = 'Betweenness';
        end
        cb.Label.String = cb_lbl;

        % Contour overlay
        hold(ax, 'on');
        valid_data = Z(isfinite(Z));
        if ~isempty(valid_data)
            vmin = min(valid_data);
            vmax = max(valid_data);
            if vmax > vmin
                step = (vmax - vmin) / 8;
                mag = 10^floor(log10(step));
                step = round(step / mag) * mag;
                step = max(step, mag);
                levels = (ceil(vmin/step)*step) : step : vmax;
                if numel(levels) > 1
                    [~, hc] = contour(ax, DV_vec, Tmax_vec, Z, levels, 'k-');
                    hc.ShowText = 'on';
                    hc.LabelSpacing = 400;
                end
            end
        end

        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);
        title(ax, ttl, 'FontSize', 11);

    % ── LCC map ──────────────────────────────────────────────────────────
    case 'lcc_map'

        lcc_size = map_data.lcc_size;   % [nDV×nTmax] integer 1-N (0=skip)
        skip_map = (lcc_size == 0);

        valid_sz = lcc_size(~skip_map);
        if isempty(valid_sz)
            min_sz = 1;  max_sz = 1;
        else
            min_sz = min(valid_sz(:));
            max_sz = max(valid_sz(:));
        end

        lcc_cmap = i_lcc_colors(N_fam);   % [N_fam×3]

        rgb = ones(nTmax, nDV, 3);
        for di = 1:nDV
            for dj = 1:nTmax
                sz = lcc_size(di, dj);
                if skip_map(di, dj) || sz < 1
                    rgb(dj, di, :) = [1 1 1];
                else
                    rgb(dj, di, :) = lcc_cmap(sz, :);
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');

        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

        hold(ax, 'on');
        n_levels = max_sz - min_sz + 1;
        hp = gobjects(n_levels, 1);
        for k = 1:n_levels
            sz_k = min_sz + k - 1;
            hp(k) = patch(ax, NaN, NaN, lcc_cmap(sz_k,:), ...
                          'EdgeColor', 'none', ...
                          'DisplayName', sprintf('%d', sz_k));
        end
        lgd = legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'on');
        lgd.Title.String = 'Largest connected component size';

        title(ax, 'Largest connected component size', 'FontSize', 11);

    % ── Articulation-point count map ─────────────────────────────────────
    case 'articulation_map'

        ap_count = double(map_data.ap_count);   % [nDV×nTmax]
        skip_map = map_data.skip_map;

        max_ap   = max(ap_count(:));
        if max_ap == 0, max_ap = 1; end

        n_colors  = max_ap + 1;
        ap_cmap   = parula(max(n_colors, 2));

        rgb = ones(nTmax, nDV, 3);
        for di = 1:nDV
            for dj = 1:nTmax
                if skip_map(di, dj)
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

        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

        hold(ax, 'on');
        n_show = min(n_colors, 8);
        hp = gobjects(n_show, 1);
        for k = 1:n_show
            ci = round((k-1) * (n_colors-1) / max(n_show-1, 1)) + 1;
            ci = max(1, min(ci, n_colors));
            hp(k) = patch(ax, NaN, NaN, ap_cmap(ci,:), 'EdgeColor', 'none', ...
                          'DisplayName', sprintf('%d', ci-1));
        end
        lgd = legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'on');
        lgd.Title.String = '# articulation pts';

        title(ax, 'Articulation-point count per snapshot', 'FontSize', 11);

    otherwise
        error('net_plot_winner_map: unknown plot_type ''%s''.', plot_type);
end

% ── Common axis labels + formatting ─────────────────────────────────────────
xlabel(ax, '\DeltaV_{cap}  [m/s]',  'FontSize', 10);
ylabel(ax, 'T_{cap}  [days]', 'FontSize', 10);
set(ax, 'FontSize', 9, 'Box', 'on', 'TickDir', 'out', 'Layer', 'top');
xlim(ax, [DV_vec(1), DV_vec(end)]);
ylim(ax, [Tmax_vec(1), Tmax_vec(end)]);

% ── Save: PDF + PNG + .fig ───────────────────────────────────────────────────
pdf_path = fullfile(out_dir, [plot_type '.pdf']);
png_path = fullfile(out_dir, [plot_type '.png']);
fig_path = fullfile(out_dir, [plot_type '.fig']);

% PDF: vector output
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8.2 6], ...
         'PaperSize', [8.2 6]);
print(fig, pdf_path, '-dpdf', '-bestfit');

% PNG: raster at 150 dpi
print(fig, png_path, '-dpng', '-r150');

% .fig: MATLAB figure file
saveas(fig, fig_path);

close(fig);

end

% ═══════════════════════════════════════════════════════════════════════════
%  Local helpers
% ═══════════════════════════════════════════════════════════════════════════

function i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map)
%I_OVERLAY_LCC_CONTOUR  Red dashed contour where lcc_full = 1 (fully connected).

hold(ax, 'on');

Z = lcc_full_map';   % [nTmax × nDV]

if any(Z(:) == 0) && any(Z(:) == 1)
    contour(ax, DV_vec, Tmax_vec, Z, [0.5 0.5], ...
            'r--', 'LineWidth', 2.0, 'DisplayName', 'LCC full');
end

end

function c = i_fam_colors(N)
%I_FAM_COLORS  N visually distinct colours for the N families.
base = lines(max(N, 7));
if N <= size(base, 1)
    c = base(1:N, :);
else
    c = hsv(N);
end
end

function c = i_lcc_colors(N)
%I_LCC_COLORS  N distinct sequential colours for LCC-size values 1..N.
c = parula(max(N, 2));
end
