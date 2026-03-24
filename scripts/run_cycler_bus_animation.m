%% RUN_CYCLER_BUS_ANIMATION  Animated GIF: Cycler 11a bus dispensing craft
%
% Renders a spacecraft 'bus' travelling on the Cycler 11a periodic orbit
% in the Earth-Moon CR3BP.  As the bus moves it dispenses smaller craft at
% each departure phase; each craft travels its 3-impulse transfer arc to a
% target orbit family.  On arrival the target PO fades in and the craft
% coasts along it briefly.  The bus keeps looping throughout the animation.
%
% Prerequisites:
%   - run_rs5_single_family_sweep.m must have been run with
%     famOrigin = 'Cycler 11a', producing a sweep_summary.mat
%   - Cached atlases for Cycler 11a and all target families must exist
%
% User knobs (edit the block below):
%   SWEEP_MAT      — full path to sweep_summary.mat (leave '' to auto-detect)
%   OUT_GIF        — output filename for the animated GIF
%   N_BUS_LOOP     — frames for one complete bus orbit  (controls bus speed)
%   N_TRANSIT      — frames for each transfer arc       (arc-length resampled)
%   N_COAST        — frames craft coasts on target PO after arrival
%   FADE_FRAMES    — frames over which target PO fades in on arrival
%   FRAME_DELAY    — seconds per GIF frame (1/fps); 0.05 → 20 fps
%   BUS_DOT_SZ     — marker size for the bus spacecraft
%   CRAFT_DOT_SZ   — marker size for dispensed craft
%   FIG_W, FIG_H   — figure size in pixels
%
% Outputs written to current directory (or the path in OUT_GIF):
%   cycler_bus_animation.gif  — looping animated GIF

clear; close all; clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));

%% ═══════════════════════ USER KNOBS ═════════════════════════════════════════

SWEEP_MAT    = '';                       % ← SET THIS if auto-detect fails
OUT_GIF      = 'cycler_bus_animation.gif';

N_BUS_LOOP   = 120;     % frames for one complete bus orbit
N_TRANSIT    = 90;      % frames per transfer arc (arc-length resampled)
N_COAST      = 35;      % frames coasting on target PO after arrival
FADE_FRAMES  = 12;      % frames for target PO fade-in
FRAME_DELAY  = 0.05;    % seconds per GIF frame  (0.05 = 20 fps)

BUS_DOT_SZ   = 14;      % bus marker size
CRAFT_DOT_SZ = 8;       % dispensed craft marker size
FIG_W        = 860;     % figure width  [px]
FIG_H        = 740;     % figure height [px]

%% ═══════════════════════ AUTO-LOCATE SWEEP_MAT ═══════════════════════════════

if isempty(SWEEP_MAT)
    root = fileparts(fileparts(mfilename('fullpath')));
    % Search two levels deep with wildcards only (a literal space in the
    % folder name "Cycler 11a" can trip up dir() on some MATLAB versions).
    % Filter by family name in MATLAB after the search.
    candidates = dir(fullfile(root, 'rs3_results', '*', '*', 'sweep_summary.mat'));
    keep = arrayfun(@(f) contains(f.folder,'Cycler 11a') || ...
                        contains(f.folder,'Cyc11a'), candidates);
    candidates = candidates(keep);
    if isempty(candidates)
        error(['run_cycler_bus_animation: sweep_summary.mat not found.\n', ...
               'Run run_rs5_single_family_sweep.m with famOrigin = ''Cycler 11a'' first,\n', ...
               'then set SWEEP_MAT to the full path of sweep_summary.mat.']);
    end
    [~, newest] = max([candidates.datenum]);
    SWEEP_MAT   = fullfile(candidates(newest).folder, candidates(newest).name);
    fprintf('Auto-detected sweep: %s\n', SWEEP_MAT);
end

% If a relative path was supplied, anchor it to the repo root so the
% script works regardless of MATLAB's current working directory.
if ~isempty(SWEEP_MAT) && ~java.io.File(SWEEP_MAT).isAbsolute()
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    SWEEP_MAT = fullfile(repo_root, SWEEP_MAT);
end

%% ═══════════════════════ LOAD SWEEP DATA ═════════════════════════════════════

fprintf('Loading sweep data ...\n');
D          = load(SWEEP_MAT, 'results', 'famOrigin', 'famTargets', 'cfg');
results    = D.results;
famOrigin  = D.famOrigin;
famTargets = D.famTargets;  %#ok<NASGU>
cfg        = D.cfg;

nT = numel(results);

%% ═══════════════════════ FILTER SUCCESSFUL TRANSFERS ════════════════════════

ok = false(nT, 1);
for k = 1:nT
    R = results{k};
    if isempty(R), continue; end
    if isfield(R,'error') && ~isempty(R.error), continue; end
    % Prefer DC-corrected; fall back to raw
    if isfield(R,'Tc') && isfield(R.Tc,'converged') && R.Tc.converged
        ok(k) = isfield(R,'traj_dc') && ~isempty(R.traj_dc) ...
                && isfield(R.traj_dc,'x') && numel(R.traj_dc.x) > 3;
    elseif isfield(R,'traj_raw') && ~isempty(R.traj_raw) ...
           && isfield(R.traj_raw,'x') && numel(R.traj_raw.x) > 3
        ok(k) = true;
    end
end

ok_idx = find(ok);
nCraft = numel(ok_idx);
fprintf('%d / %d transfers successful — animating %d craft.\n', nCraft, nT, nCraft);
if nCraft == 0
    error('No successful transfers found in %s.', SWEEP_MAT);
end

%% ═══════════════════════ LOAD BUS ORBIT (Cycler 11a) ════════════════════════

fprintf('Loading bus orbit (%s) ...\n', famOrigin);
Sbus = rs3_prepare_or_load_family(famOrigin, cfg);

if ~isfield(Sbus,'Xpo') || isempty(Sbus.Xpo)
    error('Sbus.Xpo not found — atlas does not contain dense PO data.');
end

% Close the orbit so the bus wraps cleanly
if norm(Sbus.Xpo(end,1:2) - Sbus.Xpo(1,1:2)) > 1e-6
    Sbus.Xpo(end+1,:) = Sbus.Xpo(1,:);
end

bus_mu = Sbus.mu;
bus_CJ = Sbus.CJ;

% Resample bus orbit to exactly N_BUS_LOOP equally-spaced (arc-length) points.
% This gives the bus a smooth, constant-speed visual motion.
[bus_rx, bus_ry] = local_arc_resample(Sbus.Xpo(:,1), Sbus.Xpo(:,2), N_BUS_LOOP);

%% ═══════════════════════ BUILD CRAFT ARRAY ═══════════════════════════════════

% Fixed family palette — consistent colour per family across the whole project.
KNOWN_FAMILIES = { ...
    'Lyapunov L1',          'Lyapunov L2', ...
    'Cycler 21',            'Cycler 11a', ...
    'Cycler 11b',           'Cycler 32', ...
    'Resonant 2:1 Stable',  'Resonant 2:1 Unstable', ...
    'Resonant 3:1 Stable',  'Resonant 3:1 Unstable', ...
    'Resonant 5:2 Stable',  'Resonant 5:2 Unstable', ...
    'Distant Prograde Orbit'};
FAM_PALETTE = lines(13);

fprintf('Loading target atlases ...\n');
Stgt  = cell(nCraft, 1);
craft = struct('famB',{}, 'traj',{}, 'depart_xy',{}, 'arrive_xy',{}, ...
               'patch_xy',{}, 'dv_mps',{}, 'color',{}, ...
               'arc_x',{}, 'arc_y',{}, 'coast_x',{}, 'coast_y',{});

for ki = 1:nCraft
    k = ok_idx(ki);
    R = results{k};

    % Choose best trajectory
    if isfield(R,'Tc') && isfield(R.Tc,'converged') && R.Tc.converged ...
            && isfield(R,'traj_dc') && ~isempty(R.traj_dc)
        traj = R.traj_dc;
    else
        traj = R.traj_raw;
    end

    craft(ki).famB       = R.famB;
    craft(ki).traj       = traj;
    craft(ki).depart_xy  = [traj.depart_x,  traj.depart_y];
    craft(ki).arrive_xy  = [traj.arrive_x,  traj.arrive_y];
    craft(ki).patch_xy   = [traj.patch_x,   traj.patch_y];
    craft(ki).dv_mps     = traj.DV_total_mps;

    % Assign colour from fixed palette (fall back to hsv if name not found)
    fam_idx = find(strcmpi(KNOWN_FAMILIES, strtrim(R.famB)), 1);
    if isempty(fam_idx)
        fam_idx = mod(ki - 1, 13) + 1;
    end
    craft(ki).color = FAM_PALETTE(fam_idx, :);

    % Load target atlas
    try
        S = rs3_prepare_or_load_family(R.famB, cfg);
        if isfield(S,'Xpo') && ~isempty(S.Xpo)
            if norm(S.Xpo(end,1:2) - S.Xpo(1,1:2)) > 1e-6
                S.Xpo(end+1,:) = S.Xpo(1,:);
            end
            Stgt{ki} = S;
        end
    catch ME
        warning('Could not load atlas for %s: %s', R.famB, ME.message);
    end
end

%% ═══════════════════════ RESAMPLE TRANSFER ARCS ═════════════════════════════

for ki = 1:nCraft
    traj = craft(ki).traj;
    [rx, ry] = local_arc_resample(traj.x, traj.y, N_TRANSIT + 1);
    craft(ki).arc_x = rx;
    craft(ki).arc_y = ry;
end

%% ═══════════════════════ RESAMPLE COAST ARCS ════════════════════════════════

for ki = 1:nCraft
    if isempty(Stgt{ki})
        craft(ki).coast_x = [];
        craft(ki).coast_y = [];
        continue;
    end
    tgt_xy = Stgt{ki}.Xpo(:, 1:2);
    av_pt  = craft(ki).arrive_xy;
    [~, ci] = min(sum((tgt_xy - av_pt).^2, 2));

    % Wrap-around segment of length N_COAST+1 starting at arrival index
    n_tgt = size(tgt_xy, 1);
    idx   = mod((ci - 1 : ci - 1 + N_COAST), n_tgt) + 1;
    seg   = tgt_xy(idx, :);
    [cx, cy] = local_arc_resample(seg(:,1), seg(:,2), N_COAST + 1);
    craft(ki).coast_x = cx;
    craft(ki).coast_y = cy;
end

%% ═══════════════════════ COMPUTE DEPARTURE / ARRIVAL FRAMES ═════════════════

bus_xpo_xy = Sbus.Xpo(:, 1:2);
N_xpo      = size(bus_xpo_xy, 1);

dep_frame    = zeros(nCraft, 1);
arrive_frame = zeros(nCraft, 1);

for ki = 1:nCraft
    dp = craft(ki).depart_xy;
    [~, dep_raw]     = min(sum((bus_xpo_xy - dp).^2, 2));
    dep_alpha        = (dep_raw - 1) / N_xpo;          % fractional phase [0,1)
    dep_frame(ki)    = max(1, round(dep_alpha * N_BUS_LOOP));
    arrive_frame(ki) = dep_frame(ki) + N_TRANSIT;
end

N_TOTAL = N_BUS_LOOP + N_TRANSIT + N_COAST;   % ensures all craft finish coasting

%% ═══════════════════════ BUILD FIGURE & STATIC ELEMENTS ═════════════════════

fprintf('Building figure ...\n');
fig = figure('Color','w', 'Visible','off', ...
             'Units','pixels', 'Position', [100 100 FIG_W FIG_H]);
ax  = axes('Parent', fig);

rs3_core_plot_cislunar_background(bus_CJ, bus_mu, ax);
hold(ax, 'on');
axis(ax, 'equal');

% Cycler 11a orbit — always visible throughout the animation
plot(ax, Sbus.Xpo(:,1), Sbus.Xpo(:,2), '-', ...
     'Color', [0.50 0.50 0.50], 'LineWidth', 1.6, 'HandleVisibility','off');

% Bus dot (gold/yellow) — updated every frame
h_bus = plot(ax, NaN, NaN, 'o', ...
             'Color',           [0.10 0.10 0.10], ...
             'MarkerFaceColor', [0.98 0.85 0.10], ...
             'MarkerSize',      BUS_DOT_SZ, ...
             'LineWidth',       1.5, ...
             'HandleVisibility','off');

% Pre-create per-craft graphic handles (updated via set() for performance)
h_trail = gobjects(nCraft, 1);   % growing arc trail
h_dot   = gobjects(nCraft, 1);   % craft dot in transit
h_po    = gobjects(nCraft, 1);   % target PO (fades in on arrival)
h_coast = gobjects(nCraft, 1);   % craft dot coasting on target PO
h_label = gobjects(nCraft, 1);   % family name label

for ki = 1:nCraft
    c = craft(ki).color;

    h_trail(ki) = plot(ax, NaN, NaN, '-', ...
                       'Color', c, 'LineWidth', 1.3, 'HandleVisibility','off');

    h_dot(ki) = plot(ax, NaN, NaN, 'o', ...
                     'Color',           local_lighten(c, 0.35), ...
                     'MarkerFaceColor', c, ...
                     'MarkerSize',      CRAFT_DOT_SZ, ...
                     'LineWidth',       1.0, ...
                     'HandleVisibility','off');

    % Target PO — drawn now but initially invisible
    if ~isempty(Stgt{ki})
        h_po(ki) = plot(ax, Stgt{ki}.Xpo(:,1), Stgt{ki}.Xpo(:,2), '--', ...
                        'Color', local_lighten(c, 0.8), ...   % start nearly white
                        'LineWidth', 1.2, 'Visible', 'off', ...
                        'HandleVisibility','off');
    else
        h_po(ki) = plot(ax, NaN, NaN, '--', 'Color', c, ...
                        'LineWidth', 1.2, 'HandleVisibility','off');
    end

    % Coasting dot on target PO
    h_coast(ki) = plot(ax, NaN, NaN, 'o', ...
                       'Color',           local_lighten(c, 0.35), ...
                       'MarkerFaceColor', c, ...
                       'MarkerSize',      CRAFT_DOT_SZ - 1, ...
                       'LineWidth',       1.0, ...
                       'HandleVisibility','off');

    % Family name label — centred on target PO if available, else at arrive_xy
    if ~isempty(Stgt{ki})
        lx = mean(Stgt{ki}.Xpo(:,1));
        ly = mean(Stgt{ki}.Xpo(:,2));
    else
        lx = craft(ki).arrive_xy(1);
        ly = craft(ki).arrive_xy(2);
    end
    h_label(ki) = text(ax, lx, ly, craft(ki).famB, ...
                       'Color', c, 'FontSize', 7, 'FontWeight', 'bold', ...
                       'HorizontalAlignment', 'center', ...
                       'Visible', 'off', 'HandleVisibility','off');
end

title(ax, sprintf('Cycler 11a  —  dispensing %d craft into cislunar families', nCraft), ...
      'FontSize', 10);
xlabel(ax, 'x  [nd]');
ylabel(ax, 'y  [nd]');

%% ═══════════════════════ ANIMATION LOOP ══════════════════════════════════════

fprintf('Rendering %d frames ...\n', N_TOTAL);
gif_written = false;

for fi = 1:N_TOTAL

    % ── Bus position (wraps continuously at N_BUS_LOOP period) ───────────────
    bus_idx = mod(fi - 1, N_BUS_LOOP) + 1;
    set(h_bus, 'XData', bus_rx(bus_idx), 'YData', bus_ry(bus_idx));

    % ── Update each craft ────────────────────────────────────────────────────
    for ki = 1:nCraft
        df = dep_frame(ki);
        af = arrive_frame(ki);

        if fi < df
            % ─ Not yet departed ─────────────────────────────────────────────
            set(h_trail(ki), 'XData', NaN, 'YData', NaN);
            set(h_dot(ki),   'XData', NaN, 'YData', NaN);

        elseif fi <= af
            % ─ In transit ───────────────────────────────────────────────────
            t_idx = min(fi - df + 1, numel(craft(ki).arc_x));
            set(h_trail(ki), 'XData', craft(ki).arc_x(1:t_idx), ...
                             'YData', craft(ki).arc_y(1:t_idx));
            set(h_dot(ki),   'XData', craft(ki).arc_x(t_idx), ...
                             'YData', craft(ki).arc_y(t_idx));

        else
            % ─ Arrived ──────────────────────────────────────────────────────
            % Full transit arc stays as a permanent trail
            set(h_trail(ki), 'XData', craft(ki).arc_x, ...
                             'YData', craft(ki).arc_y);
            set(h_dot(ki),   'XData', NaN, 'YData', NaN);  % transit dot gone

            coast_age = fi - af;   % 1, 2, 3, ...

            % Fade in target PO: colour lerps from near-white → full over FADE_FRAMES
            if ~isempty(Stgt{ki})
                fade_t = min(coast_age / FADE_FRAMES, 1.0);
                c_now  = local_lerp_color(craft(ki).color, fade_t);
                set(h_po(ki), 'Color', c_now, 'Visible', 'on');
            end
            set(h_label(ki), 'Visible', 'on');

            % Craft coasts on target PO for N_COAST frames, then hides
            if coast_age <= N_COAST && ~isempty(craft(ki).coast_x)
                c_idx = min(coast_age, numel(craft(ki).coast_x));
                set(h_coast(ki), 'XData', craft(ki).coast_x(c_idx), ...
                                 'YData', craft(ki).coast_y(c_idx));
            else
                set(h_coast(ki), 'XData', NaN, 'YData', NaN);
            end
        end
    end

    drawnow;

    % ── Capture frame → GIF ─────────────────────────────────────────────────
    frame         = getframe(fig);
    [im_idx, cmap] = rgb2ind(frame2im(frame), 256);
    if ~gif_written
        imwrite(im_idx, cmap, OUT_GIF, 'gif', ...
                'LoopCount', Inf, 'DelayTime', FRAME_DELAY);
        gif_written = true;
    else
        imwrite(im_idx, cmap, OUT_GIF, 'gif', ...
                'WriteMode', 'append', 'DelayTime', FRAME_DELAY);
    end

    if mod(fi, 50) == 0
        fprintf('  frame %d / %d\n', fi, N_TOTAL);
    end
end

close(fig);
fprintf('\nDone.  GIF saved to: %s\n', OUT_GIF);

%% ═══════════════════════ LOCAL HELPERS ════════════════════════════════════════

function [xr, yr] = local_arc_resample(x, y, n_out)
%LOCAL_ARC_RESAMPLE  Resample (x,y) curve to n_out points, uniform arc-length.
%  Gives visually constant dot speed regardless of physical velocity.
x  = x(:);  y = y(:);
ds = hypot(diff(x), diff(y));
s  = [0; cumsum(ds)];
% Tiny jitter prevents duplicate values that break interp1
s  = s + (0 : numel(s)-1)' * 1e-14 * max(s(end), 1);
sq = linspace(0, s(end), n_out);
xr = interp1(s, x, sq, 'pchip')';
yr = interp1(s, y, sq, 'pchip')';
end

function c2 = local_lighten(c, t)
%LOCAL_LIGHTEN  Blend colour c toward white by fraction t in [0,1].
c2 = min(1, c + t * (1 - c));
end

function c_out = local_lerp_color(c_target, t)
%LOCAL_LERP_COLOR  Lerp from near-white (t=0) to c_target (t=1).
c_white = [1 1 1];
c_out   = c_white + t * (c_target - c_white);
c_out   = max(0, min(1, c_out));
end
