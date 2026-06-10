function net_make_gif(DVmatrix_sweep, TOFmatrix_sweep, DV_cap_list, Tmax_list, ...
    fixed_dj, VU_mps, TU_days, budget_factor, short_names, out_dir, frame_delay)
%NET_MAKE_GIF  Animated GIF: edges appear as DVcap budget grows (fixed Tmax).
%
%   Circle layout of N=13 orbit families.
%   Grey lines  = edges that already existed in the previous frame.
%   Red lines   = edges newly feasible in this frame.
%   Node colour = family identity.
%
%   Sweeps DVcap from index 1 (smallest) to nDV (largest).
%
% Inputs
%   DVmatrix_sweep  {nDV×nTmax cell}  DV matrices (m/s)
%   TOFmatrix_sweep {nDV×nTmax cell}  TOF matrices (days)
%   DV_cap_list     [nDV×1]           nondimensional DV-cap values
%   Tmax_list       [nTmax×1]         nondimensional Tmax values
%   fixed_dj        scalar            column index in sweep for fixed Tmax
%   VU_mps          scalar            velocity unit (m/s)
%   TU_days         scalar            time unit (days)
%   budget_factor   scalar            physical budget multiplier (= 2)
%   short_names     {N×1}             short family labels
%   out_dir         char              output directory
%   frame_delay     scalar            seconds per GIF frame

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
N     = numel(short_names);

% Guard: fixed_dj out of range
if fixed_dj < 1 || fixed_dj > nTmax
    warning('net_make_gif: fixed_dj=%d is out of range [1,%d]; skipping GIF.', ...
            fixed_dj, nTmax);
    return
end

% ── Circle layout ────────────────────────────────────────────────────────────
theta_nodes = linspace(0, 2*pi, N+1);
theta_nodes = theta_nodes(1:N);
x_pos = cos(theta_nodes(:));
y_pos = sin(theta_nodes(:));

% Node colours
fam_colors = i_gif_fam_colors(N);

gif_path = fullfile(out_dir, sprintf('animation_Tmax%d.gif', fixed_dj));

% Pre-compute adjacency for all di at the fixed dj
A_all = false(N, N, nDV);
for di = 1:nDV
    DVcap_nd = DV_cap_list(di);
    Tmax_nd  = Tmax_list(fixed_dj);
    [A, ~, ~, ~, ~, ~, ~, skip] = net_build_graph( ...
        DVmatrix_sweep{di, fixed_dj}, TOFmatrix_sweep{di, fixed_dj}, ...
        DVcap_nd, Tmax_nd, VU_mps, TU_days, budget_factor);
    if ~skip
        A_all(:,:,di) = logical(A);
    end
end

Tmax_true = budget_factor * Tmax_list(fixed_dj) * TU_days;

fig = figure('Visible', 'off', 'Color', 'white', ...
             'Units', 'pixels', 'Position', [100 100 600 600]);

for di = 1:nDV
    clf(fig);
    ax = axes(fig);

    DVcap_true = budget_factor * DV_cap_list(di) * VU_mps;
    A_cur      = A_all(:, :, di);
    A_prev     = false(N, N);
    if di > 1, A_prev = A_all(:, :, di-1); end

    A_old = A_cur & A_prev;      % already existed before this frame
    A_new = A_cur & ~A_prev;     % newly appeared in this frame

    hold(ax, 'on');

    % Draw old edges (grey, thin)
    for i = 1:N
        for j = i+1:N
            if A_old(i,j)
                line(ax, [x_pos(i), x_pos(j)], [y_pos(i), y_pos(j)], ...
                     'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
            end
        end
    end

    % Draw new edges (red, thicker)
    for i = 1:N
        for j = i+1:N
            if A_new(i,j)
                line(ax, [x_pos(i), x_pos(j)], [y_pos(i), y_pos(j)], ...
                     'Color', [0.85 0.10 0.10], 'LineWidth', 1.8);
            end
        end
    end

    % Draw nodes
    for i = 1:N
        scatter(ax, x_pos(i), y_pos(i), 160, fam_colors(i,:), 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        r = 1.20;
        text(ax, x_pos(i)*r, y_pos(i)*r, short_names{i}, ...
             'FontSize', 7, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle');
    end

    axis(ax, 'equal', 'off');
    xlim(ax, [-1.45 1.45]);
    ylim(ax, [-1.45 1.45]);
    title(ax, sprintf('DV_{cap} = %.0f m/s   |   T_{max} = %.2f days', ...
          DVcap_true, Tmax_true), 'FontSize', 10);

    % Capture frame
    drawnow;
    frame = getframe(fig);
    img   = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if di == 1
        imwrite(imind, cm, gif_path, 'gif', ...
                'Loopcount', inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, gif_path, 'gif', ...
                'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end

close(fig);
fprintf('  GIF saved: %s\n', gif_path);

end

% ── Local helper ─────────────────────────────────────────────────────────────
function c = i_gif_fam_colors(N)
base = lines(max(N, 7));
if N <= size(base, 1)
    c = base(1:N, :);
else
    c = hsv(N);
end
end
