function tests = test_grid
%TEST_GRID  Unit tests for RS3 grid construction and binning.

tests = functiontests(localfunctions);
end

% ---- rs3_grid_make basic properties ----

function test_grid_has_required_fields(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.05;
    cfg.grid.dy = 0.05;
    cfg.grid.dtheta = deg2rad(10);
    grid3 = rs3_grid_make(cfg);
    verifyTrue(testCase, isfield(grid3, 'x_edges'));
    verifyTrue(testCase, isfield(grid3, 'y_edges'));
    verifyTrue(testCase, isfield(grid3, 'th_edges'));
    verifyTrue(testCase, isfield(grid3, 'Keep'));
    verifyTrue(testCase, isfield(grid3, 'Nx'));
    verifyTrue(testCase, isfield(grid3, 'Ny'));
    verifyTrue(testCase, isfield(grid3, 'Nth'));
end

function test_grid_symmetric_y_edges(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.05;
    cfg.grid.dy = 0.05;
    cfg.grid.dtheta = deg2rad(10);
    cfg.grid.enforce_y0_edge = true;
    grid3 = rs3_grid_make(cfg);
    % y=0 must be an edge
    verifyTrue(testCase, any(abs(grid3.y_edges) < 1e-14));
end

function test_grid_theta_full_circle(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.05;
    cfg.grid.dy = 0.05;
    cfg.grid.dtheta = deg2rad(10);
    grid3 = rs3_grid_make(cfg);
    % Theta edges should span [-pi, pi]
    verifyEqual(testCase, grid3.th_edges(1), -pi, 'AbsTol', 1e-10);
    verifyEqual(testCase, grid3.th_edges(end), pi, 'AbsTol', 1e-10);
end

% ---- rs3_bin_xyth ----

function test_bin_center_of_domain(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.1;
    cfg.grid.dy = 0.1;
    cfg.grid.dtheta = deg2rad(10);
    grid3 = rs3_grid_make(cfg);
    % Origin should be inside the grid
    [ix, iy, it] = rs3_bin_xyth(0, 0, 0, grid3);
    verifyFalse(testCase, isnan(ix));
    verifyFalse(testCase, isnan(iy));
    verifyFalse(testCase, isnan(it));
end

function test_bin_out_of_bounds(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.1;
    cfg.grid.dy = 0.1;
    cfg.grid.dtheta = deg2rad(10);
    grid3 = rs3_grid_make(cfg);
    % Way outside domain
    [ix, ~, ~] = rs3_bin_xyth(100, 100, 0, grid3);
    verifyTrue(testCase, isnan(ix));
end

function test_bin_theta_wrapping(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.1;
    cfg.grid.dy = 0.1;
    cfg.grid.dtheta = deg2rad(10);
    grid3 = rs3_grid_make(cfg);
    % Angle of 3*pi should wrap to same bin as pi (approximately)
    [~, ~, it1] = rs3_bin_xyth(0, 0, 0.5, grid3);
    [~, ~, it2] = rs3_bin_xyth(0, 0, 0.5 + 2*pi, grid3);
    verifyEqual(testCase, it1, it2);
end


function test_keep_mask_excludes_primaries(testCase)
    cfg = rs3_cfg_defaults();
    cfg.grid.dx = 0.02;
    cfg.grid.dy = 0.02;
    grid3 = rs3_grid_make(cfg);

    mu = 0.012150584270572;
    CJ = 3.1;
    Keep = rs3_keep_mask_xy(grid3, CJ, mu, cfg.sys.RE_nd, cfg.sys.RM_nd);

    % Find nearest centers to Earth and Moon centers on x-axis.
    [~, ixE] = min(abs(grid3.x_centers - (-mu)));
    [~, iy0] = min(abs(grid3.y_centers - 0));
    [~, ixM] = min(abs(grid3.x_centers - (1-mu)));

    verifyFalse(testCase, Keep(iy0, ixE));
    verifyFalse(testCase, Keep(iy0, ixM));
end
