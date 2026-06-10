function tests = test_utilities
%TEST_UTILITIES  Unit tests for RS3 utility functions.

tests = functiontests(localfunctions);
end

% ---- wrap_to_pi ----

function test_wrapToPi_in_range(testCase)
    verifyEqual(testCase, wrap_to_pi(0), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, wrap_to_pi(1), 1, 'AbsTol', 1e-14);
    verifyEqual(testCase, wrap_to_pi(-1), -1, 'AbsTol', 1e-14);
end

function test_wrapToPi_wraps_large_angles(testCase)
    verifyEqual(testCase, wrap_to_pi(3*pi), wrap_to_pi(pi), 'AbsTol', 1e-14);
    verifyEqual(testCase, wrap_to_pi(-3*pi), wrap_to_pi(-pi), 'AbsTol', 1e-14);
    verifyEqual(testCase, wrap_to_pi(2*pi), 0, 'AbsTol', 1e-14);
end

function test_wrapToPi_vector(testCase)
    a = [0 pi/2 pi 3*pi/2 2*pi];
    w = wrap_to_pi(a);
    verifySize(testCase, w, [1 5]);
    % All results must be in [-pi, pi)
    verifyGreaterThanOrEqual(testCase, w, -pi);
    verifyLessThan(testCase, w, pi);
end

% ---- circ_diff ----

function test_circ_diff_zero(testCase)
    verifyEqual(testCase, circ_diff(1, 1), 0, 'AbsTol', 1e-14);
end

function test_circ_diff_symmetric(testCase)
    d1 = circ_diff(0.1, 3.0);
    d2 = circ_diff(3.0, 0.1);
    verifyEqual(testCase, d1, -d2, 'AbsTol', 1e-14);
end

function test_circ_diff_across_boundary(testCase)
    % Going from just below pi to just above -pi should be a small step
    d = circ_diff(pi - 0.1, -pi + 0.1);
    verifyLessThan(testCase, abs(d), 0.3);
end

function test_circ_diff_output_convention(testCase)
    % Output must always be in [-pi, pi) — same convention as wrap_to_pi.
    % Specifically +pi must NOT appear in the output.
    probes = [0, pi/4, pi/2, pi, -pi, 3*pi/2, -3*pi/4, 2*pi, -2*pi];
    for a = probes
        for b = probes
            d = circ_diff(a, b);
            verifyGreaterThanOrEqual(testCase, d, -pi);
            verifyLessThan(testCase, d, pi);
        end
    end
end

function test_wrapToPi_seam(testCase)
    % pi maps to -pi (seam is at +pi, excluded)
    verifyEqual(testCase, wrap_to_pi(pi), -pi, 'AbsTol', 1e-14);
    verifyEqual(testCase, wrap_to_pi(-pi), -pi, 'AbsTol', 1e-14);
end

% ---- atlas_family_short_tag ----

function test_short_tag_lyapunov(testCase)
    verifyEqual(testCase, atlas_family_short_tag('Lyapunov L1'), 'LyapL1');
    verifyEqual(testCase, atlas_family_short_tag('Lyapunov L2'), 'LyapL2');
end

function test_short_tag_fallback_sanitizes(testCase)
    tag = atlas_family_short_tag('Cycler 21');
    verifyTrue(testCase, ~isempty(tag));
    % Should contain only alphanumeric
    verifyTrue(testCase, all(isstrprop(tag, 'alphanum')));
end

function test_short_tag_truncation(testCase)
    tag = atlas_family_short_tag('A Very Long Family Name That Exceeds Limit');
    verifyLessThanOrEqual(testCase, numel(tag), 12);
end

