function h = cr3bp_plot_background(CJ, mu, ax)
%CR3BP_PLOT_BACKGROUND  ZVC + primaries background plot.
% Reused (copied) from baseline plot_cislunar_background.
    if nargin < 3 || isempty(ax)
        figure('Color','w'); ax = axes; hold(ax,'on');
    else
        hold(ax,'on');
    end

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

    img = imagesc(ax, linspace(xlim_guess(1), xlim_guess(2), Nx), ...
        linspace(ylim_guess(1), ylim_guess(2), Ny), ones(size(Z)));
    set(img, 'AlphaData', 0.7*mask, 'AlphaDataMapping','none');
    set(img, 'CData', cat(3, 0.92*ones(size(Z)), 0.92*ones(size(Z)), 0.92*ones(size(Z))));
    set(ax, 'YDir', 'normal');
    uistack(img, 'bottom');

    contour(ax, X, Y, Z, [0 0], 'k', 'LineWidth', 1.25);

    filledCircle(ax, [xE, yE], RE, [0.2 0.45 0.85], 0.9);
    filledCircle(ax, [xM, yM], RM, [0.5 0.5 0.55], 0.9);

    plot(ax, xL1, 0, 'p', 'MarkerSize', 9, 'MarkerFaceColor',[0.85 0.2 0.2], 'MarkerEdgeColor','k', 'LineWidth',0.5);
    plot(ax, xL2, 0, 'p', 'MarkerSize', 9, 'MarkerFaceColor',[0.85 0.2 0.2], 'MarkerEdgeColor','k', 'LineWidth',0.5);

    text(xL1, 0, '  L1', 'Parent', ax, 'VerticalAlignment','bottom', 'FontSize', 10);
    text(xL2, 0, '  L2', 'Parent', ax, 'VerticalAlignment','bottom', 'FontSize', 10);

    axis(ax, 'equal'); box(ax, 'on'); grid(ax, 'on');
    xlabel(ax, 'x'); ylabel(ax, 'y');
    xlim(ax, xlim_guess); ylim(ax, ylim_guess);
    title(ax, sprintf('ZVC, Earth–Moon (\\mu = %.9f),  C_J = %.6f', mu, CJ));
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

function h = filledCircle(ax, center, radius, faceColor, faceAlpha)
    t = linspace(0, 2*pi, 200);
    xx = center(1) + radius*cos(t);
    yy = center(2) + radius*sin(t);
    h = fill(ax, xx, yy, faceColor, 'EdgeColor','k', 'LineWidth',0.6, 'FaceAlpha', faceAlpha);
end
