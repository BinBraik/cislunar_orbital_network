function net_plot_baseline(baseline_data, short_names, N_fam, out_dir)
%NET_PLOT_BASELINE  Save three separate per-family centrality bar charts.
%
%   Generates one figure each for strength, harmonic closeness, and betweenness.
%   Families are sorted from largest to smallest value in each figure.
%   Articulation-point families are marked with a star above their bar.
%
%   Saves:
%     <out_dir>/baseline_strength.{pdf,png,fig}
%     <out_dir>/baseline_harmonic_closeness.{pdf,png,fig}
%     <out_dir>/baseline_betweenness.{pdf,png,fig}
%
% Inputs
%   baseline_data  struct with fields:
%     .strength            [N×1]  strength values
%     .harmonic_closeness  [N×1]  harmonic closeness values
%     .betweenness         [N×1]  betweenness values
%     .is_articulation     [N×1 logical]
%     .DVcap_mps           scalar  (m/s)
%     .Tmax_days           scalar  (days)
%   short_names    {N×1}
%   N_fam          scalar
%   out_dir        char

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

snap_title = sprintf('\\DeltaV_{cap} = %.0f m/s,  T_{cap} = %.1f days', ...
                     baseline_data.DVcap_mps, baseline_data.Tmax_days);

metrics = { ...
    'baseline_strength',           baseline_data.strength(:),           'Strength',           [0.20 0.55 0.85]; ...
    'baseline_harmonic_closeness', baseline_data.harmonic_closeness(:), 'Harmonic closeness', [0.25 0.72 0.47]; ...
    'baseline_betweenness',        baseline_data.betweenness(:),        'Betweenness',        [0.89 0.42 0.22]; ...
    };

is_ap = baseline_data.is_articulation(:);

for m = 1:size(metrics, 1)
    fname  = metrics{m, 1};
    vals   = metrics{m, 2};
    metric_lbl = metrics{m, 3};
    bar_col    = metrics{m, 4};

    % Sort descending
    [vals_s, ord] = sort(vals, 'descend');
    names_s = short_names(ord);
    ap_s    = is_ap(ord);

    x = 1:N_fam;

    fig = figure('Visible', 'off', 'Color', 'white', ...
                 'Units', 'pixels', 'Position', [100 100 900 480]);
    ax  = axes(fig);

    bar(ax, x, vals_s, 0.68, 'FaceColor', bar_col, 'EdgeColor', 'none');
    hold(ax, 'on');

    % Mark articulation points with a star above each bar
    ap_idx = find(ap_s);
    if ~isempty(ap_idx)
        y_star = vals_s(ap_idx) + 0.04 * max(vals_s(isfinite(vals_s)));
        scatter(ax, ap_idx, y_star, 70, 'k*', 'LineWidth', 1.3, ...
                'DisplayName', 'Articulation pt');
        legend(ax, 'Location', 'northeast', 'FontSize', 8);
    end

    set(ax, 'XTick', x, 'XTickLabel', names_s, ...
            'XTickLabelRotation', 35, 'FontSize', 8.5, ...
            'Box', 'on', 'TickDir', 'out');
    ylabel(ax, metric_lbl, 'FontSize', 10);
    xlim(ax, [0.4, N_fam + 0.6]);
    grid(ax, 'on');

    title(ax, sprintf('%s — %s', metric_lbl, snap_title), ...
          'FontSize', 10, 'FontWeight', 'bold');

    base = fullfile(out_dir, fname);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 9 4.8], ...
             'PaperSize', [9 4.8]);
    print(fig, [base '.pdf'], '-dpdf', '-bestfit');
    print(fig, [base '.png'], '-dpng', '-r150');
    saveas(fig, [base '.fig']);

    close(fig);
end

end
