function atlas_grid_visual_validate(grid3, cfg, outdir)
%ATLAS_GRID_VISUAL_VALIDATE  Produce Step 2 visuals for grid + theta seam.
%
% Saves a few small figures illustrating:
% - y=0 edge location within y edges
% - (x,y) grid symmetry
% - theta seam adjacency test (points near -pi and +pi)

assert(isstruct(grid3), 'grid3 must be a struct');
assert(ischar(outdir) || isstring(outdir), 'outdir must be a path');

figvis = cfg.io.fig_visible;

% ---------------------------------------------------------------------
% Fig 1: y edges with y=0 highlighted
% ---------------------------------------------------------------------
f1 = figure('Visible', figvis);
plot(grid3.y_edges, 0*grid3.y_edges, '.-'); hold on;
yl = ylim;
plot([0 0], yl, 'k-');
title('Grid: y-edges (y=0 is an edge)');
xlabel('y edge value');
ylabel('dummy');
grid on;
io_save_figure(f1, outdir, 'step2_y_edges', cfg);
close(f1);

% ---------------------------------------------------------------------
% Fig 2: show grid lines near y=0 to visually verify no straddling
% ---------------------------------------------------------------------
f2 = figure('Visible', figvis);
xs = grid3.x_edges;
ys = grid3.y_edges;

% draw only a small window around y=0 for clarity
ywin = max(5, min(numel(ys), 25));
i0 = find(abs(ys) < 10*eps, 1, 'first');
ilo = max(1, i0 - floor(ywin/2));
ihi = min(numel(ys), i0 + floor(ywin/2));
ys2 = ys(ilo:ihi);

for k = 1:numel(xs)
    plot([xs(k) xs(k)], [ys2(1) ys2(end)], '-'); hold on;
end
for k = 1:numel(ys2)
    plot([xs(1) xs(end)], [ys2(k) ys2(k)], '-');
end
axis equal tight;
title('Grid: local (x,y) view around y=0');
xlabel('x'); ylabel('y');
grid on;
io_save_figure(f2, outdir, 'step2_xy_grid_near_y0', cfg);
close(f2);

% ---------------------------------------------------------------------
% Fig 3: theta seam bin adjacency visualization
% ---------------------------------------------------------------------
eps_th = 1e-12;
thA = -pi + eps_th;
thB =  pi - eps_th;

[~,~,itA] = grid3.bin_xyth(0,0,thA);
[~,~,itB] = grid3.bin_xyth(0,0,thB);
cA = grid3.th_centers(itA);
cB = grid3.th_centers(itB);
dc = abs(circ_diff(cA,cB));

f3 = figure('Visible', figvis);
plot(grid3.th_centers, zeros(size(grid3.th_centers)), '.'); hold on;
plot(cA, 0, 'o');
plot(cB, 0, 'o');
title(sprintf('Grid: theta seam test (dc=%.3g rad, dtheta=%.3g)', dc, grid3.dtheta));
xlabel('theta bin centers (rad)');
ylabel('dummy');
grid on;
io_save_figure(f3, outdir, 'step2_theta_seam', cfg);
close(f3);

end
