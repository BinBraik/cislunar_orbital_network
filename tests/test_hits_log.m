function tests = test_hits_log
%TEST_HITS_LOG  Unit tests for rs3_hits_log_from_traj (Phase 2 rewrite).
%
% Verifies correctness of the vectorized segment-walk: output struct shape,
% voxel-index bounds, run-length deduplication, metadata stamping, 4D-state
% conversion, and segwalk >= no-segwalk row count.

tests = functiontests(localfunctions);
end

% ---- shared setup ----

function [grid3, cfg] = local_setup()
cfg = rs3_cfg_defaults();
cfg.grid.dx     = 0.1;
cfg.grid.dy     = 0.1;
cfg.grid.dtheta = deg2rad(10);
grid3 = rs3_grid_make(cfg);
end

% ---- tests ----

function test_single_point_returns_zero_rows(testCase)
    % nPts == 1 means no segments → empty output
    [grid3, cfg] = local_setup();
    t = 0;
    X = [0, 0, 0];
    rows = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    verifyEqual(testCase, double(rows.n), 0);
end

function test_output_is_packed_struct(testCase)
    [grid3, cfg] = local_setup();
    t = [0; 0.5; 1];
    X = [0, 0, 0; 0.15, 0, 0; 0.30, 0, 0];
    rows = rs3_hits_log_from_traj(1, 2, 1, 1, t, X, grid3, cfg);
    verifyTrue(testCase, isstruct(rows));
    for f = {'n','ix','iy','it','iSeed','iHead','leg','halfFlag','t'}
        verifyTrue(testCase, isfield(rows, f{1}));
    end
    verifyGreaterThan(testCase, double(rows.n), 0);
end

function test_row_indices_in_valid_range(testCase)
    [grid3, cfg] = local_setup();
    cfg.log.segwalk.enable = false;
    t = linspace(0, 1, 20)';
    X = [linspace(-0.4, 0.4, 20)', zeros(20,1), zeros(20,1)];
    rows = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    n = double(rows.n);
    if n > 0
        verifyGreaterThanOrEqual(testCase, double(rows.ix(1:n)), 1);
        verifyLessThanOrEqual(testCase,    double(rows.ix(1:n)), grid3.Nx);
        verifyGreaterThanOrEqual(testCase, double(rows.iy(1:n)), 1);
        verifyLessThanOrEqual(testCase,    double(rows.iy(1:n)), grid3.Ny);
        verifyGreaterThanOrEqual(testCase, double(rows.it(1:n)), 1);
        verifyLessThanOrEqual(testCase,    double(rows.it(1:n)), grid3.Nth);
    end
end

function test_no_consecutive_duplicate_voxels(testCase)
    % Run-length deduplication must eliminate consecutive equal (ix,iy,it)
    [grid3, cfg] = local_setup();
    t = linspace(0, 2, 50)';
    X = [linspace(-0.3, 0.3, 50)', linspace(-0.2, 0.2, 50)', ...
         linspace(-pi/4, pi/4, 50)'];
    rows = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    n = double(rows.n);
    if n > 1
        same_ix = diff(double(rows.ix(1:n))) == 0;
        same_iy = diff(double(rows.iy(1:n))) == 0;
        same_it = diff(double(rows.it(1:n))) == 0;
        verifyFalse(testCase, any(same_ix & same_iy & same_it));
    end
end

function test_metadata_fields_stamped_correctly(testCase)
    [grid3, cfg] = local_setup();
    iSeed = 3; iHead = 7; leg = 2; halfFlag = -1;
    t = [0; 1; 2];
    X = [0, 0, 0; 0.15, 0, 0; 0.30, 0, 0];
    rows = rs3_hits_log_from_traj(iSeed, iHead, leg, halfFlag, t, X, grid3, cfg);
    n = double(rows.n);
    if n > 0
        verifyTrue(testCase, all(rows.iSeed(1:n)    == uint16(iSeed)));
        verifyTrue(testCase, all(rows.iHead(1:n)    == uint16(iHead)));
        verifyTrue(testCase, all(rows.leg(1:n)      == uint8(leg)));
        verifyTrue(testCase, all(rows.halfFlag(1:n) == int8(halfFlag)));
    end
end

function test_4d_state_same_result_as_3d(testCase)
    % Full 4-column state [x y xdot ydot] → theta = atan2(ydot,xdot).
    % A trajectory with theta=0 (xdot>0, ydot=0) should give identical rows.
    [grid3, cfg] = local_setup();
    t = [0; 1; 2];
    v = 0.3;
    X4 = [0, 0, v, 0; 0.15, 0, v, 0; 0.30, 0, v, 0];   % 4-col, theta=0
    X3 = [0, 0, 0;    0.15, 0, 0;    0.30, 0, 0];        % 3-col, theta=0
    r4 = rs3_hits_log_from_traj(1, 1, 1, 1, t, X4, grid3, cfg);
    r3 = rs3_hits_log_from_traj(1, 1, 1, 1, t, X3, grid3, cfg);
    verifyEqual(testCase, double(r4.n), double(r3.n));
    if double(r4.n) > 0
        verifyEqual(testCase, r4.ix(1:r4.n), r3.ix(1:r3.n));
        verifyEqual(testCase, r4.iy(1:r4.n), r3.iy(1:r3.n));
        verifyEqual(testCase, r4.it(1:r4.n), r3.it(1:r3.n));
    end
end

function test_segwalk_at_least_as_many_rows_as_no_segwalk(testCase)
    % With a large step, segwalk subdivides and catches more voxels
    [grid3, cfg] = local_setup();
    t = [0; 1];
    X = [0, 0, 0; 0.45, 0, 0];   % single step spanning ~4-5 voxels
    cfg.log.segwalk.enable = false;
    rows_no = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    cfg.log.segwalk.enable = true;
    rows_sw = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    verifyGreaterThanOrEqual(testCase, double(rows_sw.n), double(rows_no.n));
end

function test_out_of_bounds_trajectory_returns_zero_rows(testCase)
    % A trajectory entirely outside the grid domain should produce no rows
    [grid3, cfg] = local_setup();
    R = grid3.Rdom;
    t = [0; 1];
    X = [R*2, R*2, 0; R*3, R*3, 0];   % way outside domain
    rows = rs3_hits_log_from_traj(1, 1, 1, 1, t, X, grid3, cfg);
    verifyEqual(testCase, double(rows.n), 0);
end
