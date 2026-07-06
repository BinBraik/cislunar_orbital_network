%% RUN_CYCLER_BUS_ANIMATION_DARK
%
% "Cycler bus" explainer animation (dark theme): a mothership travels once
% around an origin periodic orbit (PO) while dispensing small spacecraft
% toward several target families, using the results of
% RUN_TRAJ_SINGLE_FAMILY_SWEEP (sweep_summary.mat).
%
% Storyboard, per target:
%   - the mothership keeps moving around the origin PO the whole time
%     (departures can overlap in flight -- concurrent, not one-at-a-time)
%   - a spacecraft departs when the mothership reaches that target's real
%     optimal departure phase, then flies the already-computed
%     departure->patch->arrival transfer (traj_dc.x/y/t), timed
%     proportionally to its real time-of-flight (rescaled to fit the video)
%   - on arrival the target PO fades in, the spacecraft coasts along it for
%     a few (video) seconds, then both fade out
%
% Requires: <ROOT_DIR>/<TAG>/<origin>/sweep_summary.mat produced by
% RUN_TRAJ_SINGLE_FAMILY_SWEEP.

clear; clc;

% ===================== USER KNOBS =====================
ROOT_DIR   = fullfile(pwd, 'rs3_results');   % <-- match cfg.io.out_root from the sweep run
TAG        = '20260318_225856';              % <-- sweep run folder (cfg.io.tag)
ORIGIN_DIR = '';                             % '' = auto-detect (must be unique under TAG)

VIDEO_LAP_SEC        = 22;    % video seconds for the mothership's one lap of the origin PO
COAST_SEC            = 2.5;   % video seconds coasting on each target PO after arrival
FADE_SEC             = 1.0;   % video seconds to fade out target PO + spacecraft
END_HOLD_SEC         = 2.0;   % trailing hold once everything has finished
COAST_FRAC_OF_PERIOD = 0.12;  % fraction of the target's own period covered while coasting

FPS           = 30;
FIG_SIZE      = [1400 1050];
VIDEO_QUALITY = 95;
% ===================== END USER KNOBS =====================

setup();

% ── Locate + load the sweep results ──────────────────────────────────────────
tag_dir = fullfile(ROOT_DIR, TAG);
if ~exist(tag_dir, 'dir')
    error('[cycler_bus] Tag folder not found:\n  %s', tag_dir);
end

if isempty(ORIGIN_DIR)
    hits = dir(fullfile(tag_dir, '*', 'sweep_summary.mat'));
    if isempty(hits)
        error('[cycler_bus] No sweep_summary.mat found under:\n  %s', tag_dir);
    elseif numel(hits) > 1
        names = strjoin({hits.folder}, '\n  ');
        error('[cycler_bus] Multiple origin folders found under %s -- set ORIGIN_DIR to disambiguate:\n  %s', tag_dir, names);
    end
    summary_path = fullfile(hits(1).folder, hits(1).name);
else
    summary_path = fullfile(tag_dir, ORIGIN_DIR, 'sweep_summary.mat');
end

fprintf('[cycler_bus] Loading %s\n', summary_path);
D = load(summary_path);
famOrigin = D.famOrigin;
results   = D.results;

valid = false(numel(results), 1);
for i = 1:numel(results)
    valid(i) = isfield(results{i}, 'traj_dc') && ~isempty(results{i}.traj_dc);
end
results = results(valid);
nT = numel(results);
if nT == 0
    error('[cycler_bus] No successful transfers in %s', summary_path);
end
fprintf('[cycler_bus] Origin: %s | %d successful target transfer(s)\n', famOrigin, nT);

% ── Bright per-target colour palette (cycles if more targets than colours) ──
PALETTE = [
    0.95 0.35 0.20;   % red-orange
    0.35 0.85 0.45;   % green
    0.30 0.55 1.00;   % blue
    1.00 0.82 0.25;   % gold
    0.75 0.35 0.95;   % purple
    0.30 0.90 0.90;   % cyan
    0.95 0.55 0.75;   % pink
    0.65 0.75 0.25;   % olive
];

BG         = [0 0 0];
TXTC       = [0.92 0.94 0.99];
MOTHER_COL = [0.95 0.96 1.00];

% ── Origin PO ─────────────────────────────────────────────────────────────
ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[muO, CJO, Tf_O, X0O] = cr3bp_family_ic(famOrigin);
N_dense_O = 6001;
solO = ode113(@(t,X) cr3bp_reduced_ode(t,X,CJO,muO,false), [0 Tf_O], X0O, ode_opts);
tO   = linspace(0, Tf_O, N_dense_O)';
XpoO = deval(solO, tO)';   % [x y th]

k_time = VIDEO_LAP_SEC / Tf_O;   % video-seconds per ND-time-unit

% ── Per-target: departure phase, colour, and target PO ──────────────────────
targets = struct('famB',{}, 'traj',{}, 'traj_t_safe',{}, 'color',{}, 'Xpo',{}, ...
    't_dep_vid',{}, 't_transfer_vid',{}, 'i_arrive',{}, 'tof_days',{});

for k = 1:nT
    R  = results{k};
    tj = R.traj_dc;

    % nearest index on the origin loop for the departure point
    dd = hypot(XpoO(:,1)-tj.depart_x, XpoO(:,2)-tj.depart_y);
    [~, i_dep] = min(dd);
    t_dep_local = tO(i_dep);

    [muB, CJB, Tf_B, X0B] = cr3bp_family_ic(R.famB);
    solB = ode113(@(t,X) cr3bp_reduced_ode(t,X,CJB,muB,false), [0 Tf_B], X0B, ode_opts);
    tB   = linspace(0, Tf_B, 2001)';
    XpoB = deval(solB, tB)';

    % nearest index on the target loop for the arrival point (coast start)
    dB = hypot(XpoB(:,1)-tj.arrive_x, XpoB(:,2)-tj.arrive_y);
    [~, i_arr] = min(dB);

    tof_nd = tj.t(end);   % ND transfer duration (depart -> arrive)

    % strictly-increasing guard for interp1 (junction between arc A and
    % reversed arc B can land on/near-duplicate time samples)
    traj_t_safe = tj.t + (0:numel(tj.t)-1)' * 1e-12 * max(tj.t(end), 1);

    targets(k).famB           = R.famB;
    targets(k).traj           = tj;
    targets(k).traj_t_safe    = traj_t_safe;
    targets(k).color          = PALETTE(mod(k-1, size(PALETTE,1)) + 1, :);
    targets(k).Xpo            = XpoB;
    targets(k).i_arrive       = i_arr;
    targets(k).t_dep_vid      = t_dep_local * k_time;
    targets(k).t_transfer_vid = tof_nd * k_time;
    targets(k).tof_days       = tj.tof_total_days;

    fprintf('  [%d/%d] -> %-28s  depart phase %5.1f%% of lap | ToF %6.1f d -> %5.1fs on screen\n', ...
        k, nT, R.famB, 100*t_dep_local/Tf_O, tj.tof_total_days, targets(k).t_transfer_vid);
end

% ── Total video duration ─────────────────────────────────────────────────────
t_end_k = zeros(nT, 1);
for k = 1:nT
    t_end_k(k) = targets(k).t_dep_vid + targets(k).t_transfer_vid + COAST_SEC + FADE_SEC;
end
VIDEO_DURATION = max([VIDEO_LAP_SEC; t_end_k]) + END_HOLD_SEC;
N_frames = round(FPS * VIDEO_DURATION);
fprintf('[cycler_bus] Video duration: %.1fs (%d frames @ %d fps)\n', VIDEO_DURATION, N_frames, FPS);

% ── Figure / background ──────────────────────────────────────────────────────
fig = figure('Color', BG, 'Visible', 'off', 'InvertHardcopy', 'off', ...
    'Units', 'pixels', 'Position', [0 0 FIG_SIZE(1) FIG_SIZE(2)]);
ax = axes('Parent', fig);
cr3bp_plot_background_dark(CJO, muO, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');
set(ax, 'GridAlpha', 0.15);

title(ax, sprintf('Cycler bus -- origin: %s', famOrigin), 'Interpreter', 'none', 'Color', TXTC);
xlabel(ax, 'x [nd]', 'Color', TXTC); ylabel(ax, 'y [nd]', 'Color', TXTC);

% Origin PO (dotted)
plot(ax, [XpoO(:,1); XpoO(1,1)], [XpoO(:,2); XpoO(1,2)], ':', ...
    'Color', [0.55 0.60 0.70], 'LineWidth', 1.4, 'HandleVisibility', 'off');

% Mothership marker
h_mother = plot(ax, NaN, NaN, 'h', 'MarkerSize', 14, 'MarkerFaceColor', MOTHER_COL, ...
    'MarkerEdgeColor', [0.15 0.15 0.18], 'LineWidth', 1.2, 'DisplayName', 'Mothership');

% Per-target handles: target PO (dashed), trail, spacecraft dot, label
h_po    = gobjects(nT, 1);
h_trail = gobjects(nT, 1);
h_dot   = gobjects(nT, 1);
h_label = gobjects(nT, 1);
for k = 1:nT
    col = targets(k).color;
    h_po(k)    = plot(ax, NaN, NaN, '--', 'Color', col, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    h_trail(k) = plot(ax, NaN, NaN, '-',  'Color', col, 'LineWidth', 1.8, 'HandleVisibility', 'off');
    h_dot(k)   = plot(ax, NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, ...
                       'MarkerEdgeColor', [0.92 0.94 0.99], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    h_label(k) = text(ax, NaN, NaN, 0, sprintf('  %s (%.0f d)', local_short_name(targets(k).famB), targets(k).tof_days), ...
        'Color', col, 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
lg = legend(ax, h_mother, {'Mothership'}, 'Location', 'northwest');
set(lg, 'TextColor', TXTC, 'Color', [0.09 0.10 0.14], 'EdgeColor', [0.40 0.44 0.52]);

% ── Video writer ──────────────────────────────────────────────────────────
out_dir  = fileparts(summary_path);
vid_path = fullfile(out_dir, sprintf('cycler_bus_%s_dark.mp4', sanitize_fname(famOrigin)));
vid = VideoWriter(vid_path, 'MPEG-4');
vid.FrameRate = FPS;
vid.Quality   = VIDEO_QUALITY;
open(vid);
fprintf('[cycler_bus] Writing %d frames to:\n  %s\n', N_frames, vid_path);

for f = 1:N_frames
    t_v = (f-1) / FPS;

    % Mothership: travels for VIDEO_LAP_SEC, then freezes at end-of-lap
    t_nd_mother = min(t_v / k_time, Tf_O);
    xy_mother = interp1(tO, XpoO(:,1:2), t_nd_mother);
    set(h_mother, 'XData', xy_mother(1), 'YData', xy_mother(2));

    for k = 1:nT
        tgt = targets(k);
        t0 = tgt.t_dep_vid;                  % dispatch
        t1 = t0 + tgt.t_transfer_vid;         % arrival
        t2 = t1 + COAST_SEC;                  % coast end
        t3 = t2 + FADE_SEC;                   % fully faded

        if t_v < t0
            % not yet dispatched
            set(h_po(k),    'XData', NaN, 'YData', NaN);
            set(h_trail(k), 'XData', NaN, 'YData', NaN);
            set(h_dot(k),   'XData', NaN, 'YData', NaN);
            set(h_label(k), 'Position', [NaN NaN 0]);

        elseif t_v < t1
            % in transit
            t_nd_local = min((t_v - t0) / k_time, tgt.traj_t_safe(end));
            xy   = interp1(tgt.traj_t_safe, [tgt.traj.x, tgt.traj.y], t_nd_local);
            mask = tgt.traj_t_safe <= t_nd_local;
            set(h_trail(k), 'XData', tgt.traj.x(mask), 'YData', tgt.traj.y(mask), 'Color', tgt.color);
            set(h_dot(k),   'XData', xy(1), 'YData', xy(2), 'MarkerFaceColor', tgt.color);
            set(h_label(k), 'Position', [xy(1) xy(2) 0], 'Color', tgt.color);

            % ease the target PO in during the final 15% of the transfer
            frac = (t_v - t0) / max(tgt.t_transfer_vid, eps);
            if frac > 0.85
                a = min(1, (frac - 0.85) / 0.15);
                set(h_po(k), 'XData', [tgt.Xpo(:,1); tgt.Xpo(1,1)], 'YData', [tgt.Xpo(:,2); tgt.Xpo(1,2)], ...
                    'Color', a*tgt.color + (1-a)*BG);
            else
                set(h_po(k), 'XData', NaN, 'YData', NaN);
            end

        elseif t_v < t2
            % coasting on the target PO
            set(h_po(k), 'XData', [tgt.Xpo(:,1); tgt.Xpo(1,1)], 'YData', [tgt.Xpo(:,2); tgt.Xpo(1,2)], 'Color', tgt.color);
            frac = (t_v - t1) / COAST_SEC;
            n_B  = size(tgt.Xpo, 1);
            step = round(frac * COAST_FRAC_OF_PERIOD * n_B);
            idx  = mod(tgt.i_arrive - 1 + (0:step), n_B) + 1;
            set(h_trail(k), 'XData', tgt.Xpo(idx,1), 'YData', tgt.Xpo(idx,2), 'Color', tgt.color);
            set(h_dot(k),   'XData', tgt.Xpo(idx(end),1), 'YData', tgt.Xpo(idx(end),2), 'MarkerFaceColor', tgt.color);
            set(h_label(k), 'Position', [tgt.Xpo(idx(end),1) tgt.Xpo(idx(end),2) 0], 'Color', tgt.color);

        elseif t_v < t3
            % fading out
            a = 1 - (t_v - t2) / FADE_SEC;
            col_fade = a*tgt.color + (1-a)*BG;
            set(h_po(k),    'Color', col_fade);
            set(h_trail(k), 'Color', col_fade);
            set(h_dot(k),   'MarkerFaceColor', col_fade, 'MarkerEdgeColor', a*[0.92 0.94 0.99] + (1-a)*BG);
            set(h_label(k), 'Color', col_fade);

        else
            % done
            set(h_po(k),    'XData', NaN, 'YData', NaN);
            set(h_trail(k), 'XData', NaN, 'YData', NaN);
            set(h_dot(k),   'XData', NaN, 'YData', NaN);
            set(h_label(k), 'Position', [NaN NaN 0]);
        end
    end

    drawnow;
    writeVideo(vid, getframe(fig));

    if mod(f, max(1, round(N_frames/20))) == 0
        fprintf('    %3d%%\n', round(100*f/N_frames));
    end
end

close(vid);
close(fig);
fprintf('[cycler_bus] Done: %s\n', vid_path);

% =============================================================================
% LOCAL HELPERS
% =============================================================================

function s = local_short_name(fam_name)
    map = { ...
        'Lyapunov L1',            'Lyap. L1';     ...
        'Lyapunov L2',            'Lyap. L2';     ...
        'Cycler 21',              'Cycler 2:1';   ...
        'Cycler 11a',             'Cycler (1,1)a';...
        'Cycler 11b',             'Cycler (1,1)b';...
        'Cycler 32',              'Cycler 3:2';   ...
        'Resonant 2to1 Stable',   'R2:1-S';       ...
        'Resonant 2to1 Unstable', 'R2:1-U';       ...
        'Resonant 3to1 Stable',   'R3:1-S';       ...
        'Resonant 3to1 Unstable', 'R3:1-U';       ...
        'Resonant 5to2 Stable',   'R5:2-S';       ...
        'Resonant 5to2 Unstable', 'R5:2-U';       ...
        'Distant Prograde Orbit', 'DPO';          ...
    };
    idx = find(strcmp(map(:,1), fam_name), 1);
    if ~isempty(idx), s = map{idx,2}; else, s = fam_name; end
end
