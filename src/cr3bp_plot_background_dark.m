function h = cr3bp_plot_background_dark(CJ, mu, ax)
%CR3BP_PLOT_BACKGROUND_DARK  Dark-theme ZVC + primaries background.
% Same geometry as cr3bp_plot_background, restyled for dark-theme
% presentation figures (poster/paper-style: near-black space, glowing
% primaries, light ZVC boundary).

if nargin < 3 || isempty(ax)
    figure('Color', [0.047 0.055 0.078]); ax = axes; hold(ax,'on');
else
    hold(ax,'on');
end

BG      = [0.047 0.055 0.078];   % near-black navy "space"
ALLOWED = [0.11  0.15  0.24];    % Hill-allowed region tint
ZVCLINE = [0.60  0.88  1.00];    % light cyan ZVC boundary
EARTHC  = [0.30  0.58  1.00];
MOONC   = [0.78  0.78  0.83];
LPTC    = [1.00  0.82  0.25];    % gold Lagrange markers
TXTC    = [0.90  0.93  0.98];
GRIDC   = [0.32  0.36  0.45];

xE = -mu; yE = 0;
xM = 1 - mu; yM = 0;

RE = 6378.0/384400.0;
RM = 1737.0/384400.0;

Ux_on_xaxis = @(x) x ...
    - (1-mu)*(x+mu)/abs(x+mu)^3 ...
    - mu*(x-1+mu)/abs(x-1+mu)^3;

gamma = (mu/3)^(1/3);
xL1_0 = xM - gamma;
xL2_0 = xM + gamma;

opts = optimset('TolX',1e-14,'Display','off');
xL1 = fzero(Ux_on_xaxis, xL1_0, opts);
xL2 = fzero(Ux_on_xaxis, xL2_0, opts);

xlim_guess = [-1.4, 1.4];
ylim_guess = [-1.4, 1.4];

Nx = 800; Ny = 600;
[X, Y] = meshgrid(linspace(xlim_guess(1), xlim_guess(2), Nx), ...
                  linspace(ylim_guess(1), ylim_guess(2), Ny));

[U, ~, ~] = U_and_grad(X, Y, mu);
Z = 2*U - CJ;
mask = Z < 0;

set(ax, 'Color', BG);

img = imagesc(ax, linspace(xlim_guess(1), xlim_guess(2), Nx), ...
    linspace(ylim_guess(1), ylim_guess(2), Ny), ones(size(Z)));
set(img, 'AlphaData', 0.85*mask, 'AlphaDataMapping','none');
set(img, 'CData', cat(3, ALLOWED(1)*ones(size(Z)), ALLOWED(2)*ones(size(Z)), ALLOWED(3)*ones(size(Z))));
set(ax, 'YDir', 'normal');
uistack(img, 'bottom');

contour(ax, X, Y, Z, [0 0], 'Color', ZVCLINE, 'LineWidth', 1.35);

filledGlowCircle(ax, [xE, yE], RE, EARTHC);
filledGlowCircle(ax, [xM, yM], RM, MOONC);

plot(ax, xL1, 0, 'p', 'MarkerSize', 9, 'MarkerFaceColor', LPTC, 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);
plot(ax, xL2, 0, 'p', 'MarkerSize', 9, 'MarkerFaceColor', LPTC, 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);

text(xL1, 0, '  L1', 'Parent', ax, 'VerticalAlignment','bottom', 'FontSize', 10, 'Color', TXTC);
text(xL2, 0, '  L2', 'Parent', ax, 'VerticalAlignment','bottom', 'FontSize', 10, 'Color', TXTC);

axis(ax, 'equal'); box(ax, 'on'); grid(ax, 'on');
set(ax, 'GridColor', GRIDC, 'GridAlpha', 0.4, 'XColor', TXTC, 'YColor', TXTC);
xlabel(ax, 'x', 'Color', TXTC); ylabel(ax, 'y', 'Color', TXTC);
xlim(ax, xlim_guess); ylim(ax, ylim_guess);
title(ax, sprintf('ZVC, Earth-Moon (\\mu = %.9f),  C_J = %.6f', mu, CJ), 'Color', TXTC);
h = struct();
end

function [U, Ux, Uy] = U_and_grad(x, y, mu)
    x1 = x + mu;    y1 = y;
    x2 = x - 1 + mu; y2 = y;
    r1 = sqrt(x1.^2 + y1.^2);
    r2 = sqrt(x2.^2 + y2.^2);
    U  = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2;
    Ux = x - (1-mu).*x1./(r1.^3) - mu.*x2./(r2.^3);
    Uy = y - (1-mu).*y1./(r1.^3) - mu.*y2./(r2.^3);
end

function h = filledGlowCircle(ax, center, radius, faceColor)
% Solid disk plus a couple of soft, semi-transparent halo rings for a
% "glow" look against the dark background.
    t = linspace(0, 2*pi, 200);
    for k = 3:-1:1
        rr = radius * (1 + 0.9*k);
        aa = 0.10 / k;
        fill(ax, center(1)+rr*cos(t), center(2)+rr*sin(t), faceColor, ...
            'EdgeColor', 'none', 'FaceAlpha', aa, 'HandleVisibility', 'off');
    end
    h = fill(ax, center(1)+radius*cos(t), center(2)+radius*sin(t), faceColor, ...
        'EdgeColor', [0.95 0.97 1.00], 'LineWidth', 0.6, 'FaceAlpha', 0.95);
end
