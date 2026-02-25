function tests = test_rs4_dc_jacobian
%TEST_RS4_DC_JACOBIAN  Unit tests for RS4 DC Jacobian assembly helper.

tests = functiontests(localfunctions);
end

function test_mode_auto_uses_stm_or_fd(testCase)
cfg = rs3_cfg_defaults();
cfg.rs4.dc.jacobian_mode = 'auto';

[mu, CJ, ~, X0] = rs3_core_family_ic('Lyapunov L1');
arcA = struct('IC', X0, 'tfinal', 0.05, 'CJ', CJ, 'mu', mu, 'mirror_heading', false);
arcB = struct('IC', X0, 'tfinal', 0.07, 'CJ', CJ, 'mu', mu, 'mirror_heading', true);

J = rs4_dc_assemble_heading_jacobian(arcA, arcB, cfg);
verifyTrue(testCase, ismember(J.mode_used, {'stm','fd'}));
verifySize(testCase, J.jacobian, [1 2]);
verifyTrue(testCase, all(isfinite(J.jacobian)));
verifyTrue(testCase, isfinite(J.residual_heading_rad));
end

function test_mode_fd_forces_fd(testCase)
cfg = rs3_cfg_defaults();
cfg.rs4.dc.jacobian_mode = 'fd';

[mu, CJ, ~, X0] = rs3_core_family_ic('Lyapunov L1');
arcA = struct('IC', X0, 'tfinal', 0.03, 'CJ', CJ, 'mu', mu, 'mirror_heading', false);
arcB = struct('IC', X0, 'tfinal', 0.04, 'CJ', CJ, 'mu', mu, 'mirror_heading', true);

J = rs4_dc_assemble_heading_jacobian(arcA, arcB, cfg);
verifyEqual(testCase, J.mode_used, 'fd');
verifyTrue(testCase, all(isfinite(J.jacobian)));
end

function test_mode_stm_errors_when_unavailable(testCase)
cfg = rs3_cfg_defaults();
cfg.rs4.dc.jacobian_mode = 'stm';

arcA = struct('IC', [1;0;0], 'tfinal', -1.0, 'CJ', 3.1, 'mu', 0.01215, 'mirror_heading', false);
arcB = struct('IC', [1;0;0], 'tfinal', 0.1,  'CJ', 3.1, 'mu', 0.01215, 'mirror_heading', false);

verifyError(testCase, @() rs4_dc_assemble_heading_jacobian(arcA, arcB, cfg), 'rs4:dc:stmFailed');
end


function test_mode_auto_falls_back_to_fd_when_stm_threshold_exceeded(testCase)
cfg = rs3_cfg_defaults();
cfg.rs4.dc.jacobian_mode = 'auto';
cfg.rs4.dc.stm_max_sensitivity = 1e-12;

[mu, CJ, ~, X0] = rs3_core_family_ic('Lyapunov L1');
arcA = struct('IC', X0, 'tfinal', 0.05, 'CJ', CJ, 'mu', mu, 'mirror_heading', false);
arcB = struct('IC', X0, 'tfinal', 0.07, 'CJ', CJ, 'mu', mu, 'mirror_heading', true);

J = rs4_dc_assemble_heading_jacobian(arcA, arcB, cfg);
verifyEqual(testCase, J.mode_used, 'fd');
verifyTrue(testCase, isfield(J.details, 'stm'));
verifyTrue(testCase, contains(lower(J.details.stm.reason), 'exceeded threshold'));
end
