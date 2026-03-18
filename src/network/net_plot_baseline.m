function net_plot_baseline(baseline_data, short_names, N_fam, out_dir)
%NET_PLOT_BASELINE  Bar chart of per-family centrality values at one snapshot.
%
%   Saves <out_dir>/baseline_metrics.{pdf,png,fig}.
%
% Inputs
%   baseline_data  struct with fields:
%     .strength            [N×1]  strength values
%     .harmonic_closeness  [N×1]  harmonic closeness values
%     .betweenness         [N×1]  betweenness values
%     .R_budget            [N×1]  budget reachability fraction (0-1)
%     .is_articulation     [N×1 logical]  articulation-point flag
%     .DVcap_mps           scalar  DV budget used (m/s)
%     .Tmax_days           scalar  Tmax budget used (days)
%   short_names    {N×1}  short family labels
%   N_fam          scalar  number of families
%   out_dir        char    output directory (created if absent)

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

str_vals = baseline_data.strength(:);
hc_vals  = baseline_data.harmonic_closeness(:);
bw_vals  = baseline_data.betweenness(:);
rb_vals  = baseline_data.R_budget(:);
is_ap    = baseline_data.is_articulation(:);
DVcap    = baseline_data.DVcap_mps;
Tmax     = baseline_data.Tmax_days;

% Sort families by strength (descending) for readability
[~, ord] = sort(str_vals, 'descend');

names_sorted = short_names(ord);
str_s  = str_vals(ord);
hc_s   = hc_vals(ord);
bw_s   = bw_vals(ord);
rb_s   = rb_vals(ord);
ap_s   = is_ap(ord);

x = 1:N_fam;

% ── Figure layout: 3 metric subplots + 1 reachability subplot ────────────────
fig = figure('Visible', 'off', 'Color', 'white', ...
             'Units', 'pixels', 'Position', [100 100 1100 700]);

metrics = { ...
    str_s,  'Strength',            [0.20 0.55 0.85]; ...
    hc_s,   'Harmonic closeness',  [0.25 0.72 0.47]; ...
    bw_s,   'Betweenness',         [0.89 0.42 0.22]; ...
    rb_s,   'Budget reachability', [0.62 0.42 0.78]; ...
    };

n_sub = size(metrics, 1);

for m = 1:n_sub
    ax = subplot(2, 2, m, 'Parent', fig);

    vals  = metrics{m, 1};
    lbl   = metrics{m, 2};
    col   = metrics{m, 3};

    b = bar(ax, x, vals, 0.65, 'FaceColor', col, 'EdgeColor', 'none');

    % Mark articulation points with a star marker
    ap_idx = find(ap_s);
    if ~isempty(ap_idx)
        hold(ax, 'on');
        scatter(ax, ap_idx, vals(ap_idx) * 1.05, 60, 'k*', ...
                'LineWidth', 1.2, 'DisplayName', 'Articulation pt');
    end

    set(ax, 'XTick', x, 'XTickLabel', names_sorted, ...
            'XTickLabelRotation', 35, 'FontSize', 7.5, ...
            'Box', 'on', 'TickDir', 'out');
    ylabel(ax, lbl, 'FontSize', 9);
    xlim(ax, [0.4, N_fam + 0.6]);
    grid(ax, 'on');
    grid(ax, 'minor');

    if m == n_sub
        legend(ax, 'Location', 'northeast', 'FontSize', 7);
    end
end

% Main title
sgtitle(fig, sprintf('Network centrality at \\DeltaV_{cap} = %.0f m/s,  T_{cap} = %.1f days', ...
        DVcap, Tmax), 'FontSize', 11, 'FontWeight', 'bold');

% ── Save: PDF + PNG + .fig ───────────────────────────────────────────────────
base_name = fullfile(out_dir, 'baseline_metrics');

set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 11 7], ...
         'PaperSize', [11 7]);
print(fig, [base_name '.pdf'], '-dpdf', '-bestfit');
print(fig, [base_name '.png'], '-dpng', '-r150');
saveas(fig, [base_name '.fig']);

close(fig);

end
