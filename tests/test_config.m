function tests = test_config
%TEST_CONFIG  Unit tests for RS3 configuration system.

tests = functiontests(localfunctions);
end

% ---- rs3_cfg_defaults ----

function test_defaults_returns_struct(testCase)
    cfg = rs3_cfg_defaults();
    verifyTrue(testCase, isstruct(cfg));
end

function test_defaults_has_all_sections(testCase)
    cfg = rs3_cfg_defaults();
    verifyTrue(testCase, isfield(cfg, 'sys'));
    verifyTrue(testCase, isfield(cfg, 'units'));
    verifyTrue(testCase, isfield(cfg, 'families'));
    verifyTrue(testCase, isfield(cfg, 'seed'));
    verifyTrue(testCase, isfield(cfg, 'grid'));
    verifyTrue(testCase, isfield(cfg, 'fan'));
    verifyTrue(testCase, isfield(cfg, 'propag'));
    verifyTrue(testCase, isfield(cfg, 'log'));
    verifyTrue(testCase, isfield(cfg, 'cache'));
    verifyTrue(testCase, isfield(cfg, 'overlap'));
    verifyTrue(testCase, isfield(cfg, 'rs4'));
    verifyTrue(testCase, isfield(cfg.rs4, 'dc'));
    verifyTrue(testCase, isfield(cfg, 'par'));
    verifyTrue(testCase, isfield(cfg, 'io'));
    verifyTrue(testCase, isfield(cfg, 'diag'));
end

function test_defaults_validate_passes(testCase)
    cfg = rs3_cfg_defaults();
    % Should not throw
    rs3_cfg_validate(cfg);
    verifyTrue(testCase, true);
end

function test_defaults_units_consistent(testCase)
    cfg = rs3_cfg_defaults();
    % VU = LU / TU
    VU_expected = cfg.units.LU_m / cfg.units.TU_s;
    verifyEqual(testCase, cfg.units.VU_mps, VU_expected, 'RelTol', 1e-10);
end

% ---- rs3_core_family_ic ----

function test_family_ic_all_families(testCase)
    families = {'Lyapunov L1', 'Lyapunov L2', 'Cycler 21', 'Cycler 11a', ...
                'Cycler 11b', 'Cycler 32', 'Resonant 2to1 Stable', ...
                'Resonant 2to1 Unstable', 'Resonant 3to1 Stable', ...
                'Resonant 3to1 Unstable', 'Resonant 5to2 Stable', ...
                'Resonant 5to2 Unstable', 'Distant Prograde Orbit'};
    for i = 1:numel(families)
        [mu, CJ, Tf, X0] = rs3_core_family_ic(families{i});
        verifyGreaterThan(testCase, mu, 0);
        verifyGreaterThan(testCase, CJ, 0);
        verifyGreaterThan(testCase, Tf, 0);
        verifySize(testCase, X0, [3 1]);
    end
end

function test_family_ic_unknown_errors(testCase)
    threw = false;
    try
        rs3_core_family_ic('Nonexistent');
    catch
        threw = true;
    end
    verifyTrue(testCase, threw);
end


function test_dc_jacobian_mode_validation(testCase)
    cfg = rs3_cfg_defaults();
    cfg.rs4.dc.jacobian_mode = 'bad_mode';
    verifyError(testCase, @() rs3_cfg_validate(cfg), 'rs3:cfg:badValue');
end
