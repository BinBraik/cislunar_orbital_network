function atlas_family_visual_validate(S, cfg, outdir)
%RS3_FAMILY_VISUAL_VALIDATE  Step 3 visuals: background + PO + seeds + Keep mask.
%
% Saves PNG figures into outdir.

safeName = regexprep(S.name,'[^A-Za-z0-9]+','_');

% --- Figure 1: PO and seeds on ZVC background ---
f1 = figure('Visible', cfg.io.fig_visible, 'Color','w');
ax1 = axes('Parent', f1); hold(ax1,'on');
cr3bp_plot_background(S.CJ, S.mu, ax1);
% PO trace is optional (it may be stripped from cache for memory reasons)
hasPO = isfield(S,'Xpo') && ~isempty(S.Xpo);
if hasPO
    plot(ax1, S.Xpo(:,1), S.Xpo(:,2), '-', 'LineWidth', 1.1);
end
if ~isempty(S.SeedsUpper)
    plot(ax1, S.SeedsUpper(:,1), S.SeedsUpper(:,2), '.', 'MarkerSize', 9);
end
if ~isempty(S.SeedsLower)
    plot(ax1, S.SeedsLower(:,1), S.SeedsLower(:,2), '.', 'MarkerSize', 6);
end
if hasPO
    title(ax1, sprintf('rs3 Step3: %s | PO + seeds', S.name), 'Interpreter','none');
else
    title(ax1, sprintf('rs3 Step3: %s | seeds (PO not cached)', S.name), 'Interpreter','none');
end
drawnow;
saveas(f1, fullfile(outdir, sprintf('step3_%s_po_seeds.png', safeName)));
close(f1);

% --- Figure 2: Keep mask preview (xy) ---
f2 = figure('Visible', cfg.io.fig_visible, 'Color','w');
ax2 = axes('Parent', f2);
imagesc(ax2, S.grid3.x_centers, S.grid3.y_centers, double(S.grid3.Keep));
set(ax2, 'YDir','normal'); axis(ax2,'equal'); grid(ax2,'on');
xlabel(ax2,'x'); ylabel(ax2,'y');
title(ax2, sprintf('rs3 Step3: %s | Keep mask (CJ allowed)', S.name), 'Interpreter','none');
hold(ax2,'on');
if hasPO
    plot(ax2, S.Xpo(:,1), S.Xpo(:,2), 'k-', 'LineWidth', 0.8);
    plot(ax2, S.Xpo(:,1), -S.Xpo(:,2), 'k-', 'LineWidth', 0.8);
end
drawnow;
saveas(f2, fullfile(outdir, sprintf('step3_%s_keep_mask.png', safeName)));
close(f2);

% --- Figure 3: theta distribution for seeds (upper) ---
f3 = figure('Visible', cfg.io.fig_visible, 'Color','w');
ax3 = axes('Parent', f3);
if isempty(S.SeedsUpper)
    text(0.1,0.5,'No seeds','Parent',ax3);
    axis(ax3,'off');
else
    histogram(ax3, S.SeedsUpper(:,3), 50);
    xlabel(ax3,'theta (rad)'); ylabel(ax3,'count');
    title(ax3, sprintf('rs3 Step3: %s | upper seed theta', S.name), 'Interpreter','none');
    grid(ax3,'on');
end
saveas(f3, fullfile(outdir, sprintf('step3_%s_theta_hist.png', safeName)));
close(f3);

end
