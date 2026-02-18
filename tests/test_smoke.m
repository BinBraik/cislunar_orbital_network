function tests = test_smoke
%TEST_SMOKE  Minimal rs3 smoke tests (fast, no heavy integrations).
%
% Run with:
%   rs3_setup
%   results = runtests('tests');
%   assertSuccess(results)

tests = functiontests(localfunctions);
end

function test_cfg_defaults_validate(testCase)
cfg = rs3_cfg_defaults();
rs3_cfg_validate(cfg);
verifyTrue(testCase, isstruct(cfg));
end

function test_family_ic_contract(testCase)
[mu,CJ,Tf,X0] = rs3_core_family_ic('Lyapunov L1');
verifyGreaterThan(testCase, mu, 0);
verifyGreaterThan(testCase, Tf, 0);
verifyGreaterThan(testCase, CJ, 0);
verifySize(testCase, X0, [3 1]);
end

function test_grid_make_validate(testCase)
cfg = rs3_cfg_defaults();
% coarsen for speed
cfg.grid.dx = 0.05;
cfg.grid.dy = 0.05;
cfg.grid.dtheta = deg2rad(5);

grid3 = rs3_grid_make(cfg);
rs3_grid_validate(grid3, cfg);
verifyTrue(testCase, isfield(grid3,'Keep'));
end
