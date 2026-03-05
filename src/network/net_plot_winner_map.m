function net_plot_winner_map(plot_type, map_data, DV_vec, Tmax_vec, ...
    lcc_full_map, short_names, N_fam, out_dir)
%NET_PLOT_WINNER_MAP  Create, save (PNG + .fig) one network-centrality figure.
%
%   Called once per figure type.  Saves <out_dir>/<plot_type>.{png,fig}.
%
% Inputs
%   plot_type     char    one of: 'winner_map_reach' | 'winner_map_gateway'
%                         | 'winner_map_hub' | 'lcc_map' | 'articulation_map'
%   map_data      struct  fields depend on plot_type (see below)
%   DV_vec        [nDV×1]   x-axis values (DVcap_true_mps, m/s)
%   Tmax_vec      [nTmax×1] y-axis values (Tmax_true_days, days)
%   lcc_full_map  [nDV×nTmax] binary; 1 where entire graph is connected (N nodes)
%   short_names   {N_fam×1} short family labels
%   N_fam         scalar    number of families (13)
%   out_dir       char      output directory (created if absent)
%
% map_data fields by plot_type
%   winner_map_*  .winner_idx  [nDV×nTmax] integer 1-N, 0=Tie, NaN=skip
%                 .tie_size    [nDV×nTmax]
%   lcc_map       .lcc_size    [nDV×nTmax] integer 1-N, 0 if skip
%   articulation_map .ap_count [nDV×nTmax] integer 0-(N-1)
%                    .skip_map [nDV×nTmax] logical

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

nDV   = numel(DV_vec);
nTmax = numel(Tmax_vec);

% ── Build figure ─────────────────────────────────────────────────────────────
fig = figure('Visible', 'off', 'Color', 'white', ...
             'Units', 'pixels', 'Position', [100 100 820 600]);
ax  = axes(fig);

switch plot_type

    % ── Winner maps (reach / gateway / hub) ──────────────────────────────
    case {'winner_map_reach', 'winner_map_gateway', 'winner_map_hub'}

        winner_idx = map_data.winner_idx;   % [nDV×nTmax]

        % 13 distinct family colours (lines(13) palette)
        fam_colors = i_fam_colors(N_fam);
        tie_color  = [0.65 0.65 0.65];      % grey

        % Build RGB image [nTmax × nDV × 3]  (imagesc row=y, col=x)
        rgb = ones(nTmax, nDV, 3);           % default white = skip

        for di = 1:nDV
            for dj = 1:nTmax
                idx = winner_idx(di, dj);
                if isnan(idx)
                    rgb(dj, di, :) = [1 1 1];              % skip → white
                elseif idx == 0
                    rgb(dj, di, :) = tie_color;            % tie  → grey
                else
                    rgb(dj, di, :) = fam_colors(idx, :);  % family colour
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');

        % Overlay lcc_full contour
        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

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

        legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, ...
               'Box', 'on');

        % Metric-specific title
        switch plot_type
            case 'winner_map_reach'
                ttl = 'Best Launch Hub (REACH winner)';
            case 'winner_map_gateway'
                ttl = 'Critical Relay Node (GATEWAY winner)';
            case 'winner_map_hub'
                ttl = 'Best Direct-Transfer Node (HUB winner)';
        end
        title(ax, ttl, 'FontSize', 11);

    % ── LCC map ──────────────────────────────────────────────────────────
    case 'lcc_map'

        lcc_size = map_data.lcc_size;   % [nDV×nTmax] integer 1-N (0=skip)
        skip_map = (lcc_size == 0);

        % N distinct colours for integer values 1..N_fam
        lcc_cmap = i_lcc_colors(N_fam);   % [N_fam×3]

        % Build RGB image
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

        % Overlay lcc_full contour (on the LCC map itself too)
        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

        % Colourbar: patches with integer labels 1..N_fam
        hold(ax, 'on');
        hp = gobjects(N_fam, 1);
        for k = 1:N_fam
            hp(k) = patch(ax, NaN, NaN, lcc_cmap(k,:), ...
                          'EdgeColor', 'none', ...
                          'DisplayName', sprintf('%d', k));
        end
        legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, ...
               'Title', 'LCC size (# nodes)', 'Box', 'on');

        title(ax, 'Largest Connected Component Size', 'FontSize', 11);

    % ── Articulation-point count map ─────────────────────────────────────
    case 'articulation_map'

        ap_count = double(map_data.ap_count);   % [nDV×nTmax]
        skip_map = map_data.skip_map;

        max_ap   = max(ap_count(:));
        if max_ap == 0, max_ap = 1; end  % avoid degenerate colourmap

        % Perceptually uniform colourmap (parula) for 0..max_ap
        n_colors  = max_ap + 1;          % 0, 1, …, max_ap
        ap_cmap   = parula(max(n_colors, 2));

        % Build RGB image
        rgb = ones(nTmax, nDV, 3);
        for di = 1:nDV
            for dj = 1:nTmax
                if skip_map(di, dj)
                    rgb(dj, di, :) = [1 1 1];
                else
                    ci = ap_count(di, dj) + 1;  % 1-based index
                    ci = max(1, min(ci, n_colors));
                    rgb(dj, di, :) = ap_cmap(ci, :);
                end
            end
        end

        image(ax, DV_vec, Tmax_vec, rgb);
        axis(ax, 'xy');

        % Overlay lcc_full contour
        i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map);

        % Legend
        hold(ax, 'on');
        n_show = min(n_colors, 8);    % cap legend entries for readability
        hp = gobjects(n_show, 1);
        for k = 1:n_show
            ci = round((k-1) * (n_colors-1) / max(n_show-1, 1)) + 1;
            ci = max(1, min(ci, n_colors));
            hp(k) = patch(ax, NaN, NaN, ap_cmap(ci,:), 'EdgeColor', 'none', ...
                          'DisplayName', sprintf('%d', ci-1));
        end
        legend(ax, hp, 'Location', 'eastoutside', 'FontSize', 7, ...
               'Title', '# articulation pts', 'Box', 'on');

        title(ax, 'Articulation-Point Count per Snapshot', 'FontSize', 11);

    otherwise
        error('net_plot_winner_map: unknown plot_type ''%s''.', plot_type);
end

% ── Common axis labels + formatting ─────────────────────────────────────────
xlabel(ax, 'DV_{cap} budget  [m/s]',  'FontSize', 10);
ylabel(ax, 'T_{max} budget  [days]',  'FontSize', 10);
set(ax, 'FontSize', 9, 'Box', 'on', 'TickDir', 'out', 'Layer', 'top');
xlim(ax, [DV_vec(1), DV_vec(end)]);
ylim(ax, [Tmax_vec(1), Tmax_vec(end)]);

% ── Save ─────────────────────────────────────────────────────────────────────
png_path = fullfile(out_dir, [plot_type '.png']);
fig_path = fullfile(out_dir, [plot_type '.fig']);

print(fig, png_path, '-dpng', '-r150');
saveas(fig, fig_path);

close(fig);

end

% ═══════════════════════════════════════════════════════════════════════════
%  Local helpers
% ═══════════════════════════════════════════════════════════════════════════

function i_overlay_lcc_contour(ax, DV_vec, Tmax_vec, lcc_full_map)
%I_OVERLAY_LCC_CONTOUR  Bold black contour where lcc_full = 1 on axes ax.

hold(ax, 'on');

% contour(x, y, Z) expects Z as [numel(y) × numel(x)]
Z = lcc_full_map';   % [nTmax × nDV]

% Only draw contour if the map has both 0 and 1 values
if any(Z(:) == 0) && any(Z(:) == 1)
    contour(ax, DV_vec, Tmax_vec, Z, [0.5 0.5], ...
            'k-', 'LineWidth', 2.5, 'DisplayName', 'LCC full');
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
%I_LCC_COLORS  N distinct colours for LCC-size values 1..N.
% Use a fixed, perceptually varied palette.
c = lines(max(N, 7));
if N > size(c, 1)
    c = hsv(N);
else
    c = c(1:N, :);
end
end
