function tests = test_utilities
%TEST_UTILITIES  Unit tests for RS3 utility functions.

tests = functiontests(localfunctions);
end

% ---- rs3_wrapToPi ----

function test_wrapToPi_in_range(testCase)
    verifyEqual(testCase, rs3_wrapToPi(0), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, rs3_wrapToPi(1), 1, 'AbsTol', 1e-14);
    verifyEqual(testCase, rs3_wrapToPi(-1), -1, 'AbsTol', 1e-14);
end

function test_wrapToPi_wraps_large_angles(testCase)
    verifyEqual(testCase, rs3_wrapToPi(3*pi), rs3_wrapToPi(pi), 'AbsTol', 1e-14);
    verifyEqual(testCase, rs3_wrapToPi(-3*pi), rs3_wrapToPi(-pi), 'AbsTol', 1e-14);
    verifyEqual(testCase, rs3_wrapToPi(2*pi), 0, 'AbsTol', 1e-14);
end

function test_wrapToPi_vector(testCase)
    a = [0 pi/2 pi 3*pi/2 2*pi];
    w = rs3_wrapToPi(a);
    verifySize(testCase, w, [1 5]);
    % All results must be in [-pi, pi)
    verifyGreaterThanOrEqual(testCase, w, -pi);
    verifyLessThan(testCase, w, pi);
end

% ---- rs3_circ_diff ----

function test_circ_diff_zero(testCase)
    verifyEqual(testCase, rs3_circ_diff(1, 1), 0, 'AbsTol', 1e-14);
end

function test_circ_diff_symmetric(testCase)
    d1 = rs3_circ_diff(0.1, 3.0);
    d2 = rs3_circ_diff(3.0, 0.1);
    verifyEqual(testCase, d1, -d2, 'AbsTol', 1e-14);
end

function test_circ_diff_across_boundary(testCase)
    % Going from just below pi to just above -pi should be a small step
    d = rs3_circ_diff(pi - 0.1, -pi + 0.1);
    verifyLessThan(testCase, abs(d), 0.3);
end

% ---- rs3_family_short_tag ----

function test_short_tag_lyapunov(testCase)
    verifyEqual(testCase, rs3_family_short_tag('Lyapunov L1'), 'LyapL1');
    verifyEqual(testCase, rs3_family_short_tag('Lyapunov L2'), 'LyapL2');
end

function test_short_tag_fallback_sanitizes(testCase)
    tag = rs3_family_short_tag('Cycler 21');
    verifyTrue(testCase, ~isempty(tag));
    % Should contain only alphanumeric
    verifyTrue(testCase, all(isstrprop(tag, 'alphanum')));
end

function test_short_tag_truncation(testCase)
    tag = rs3_family_short_tag('A Very Long Family Name That Exceeds Limit');
    verifyLessThanOrEqual(testCase, numel(tag), 12);
end

% ---- rs3_group_sorted_ids ----

function test_group_sorted_simple(testCase)
    [u, s, e] = rs3_group_sorted_ids([1; 1; 2; 2; 2; 3]);
    verifyEqual(testCase, u, [1; 2; 3]);
    verifyEqual(testCase, s, [1; 3; 6]);
    verifyEqual(testCase, e, [2; 5; 6]);
end

function test_group_sorted_empty(testCase)
    [u, s, e] = rs3_group_sorted_ids([]);
    verifyEmpty(testCase, u);
    verifyEmpty(testCase, s);
    verifyEmpty(testCase, e);
end

function test_group_sorted_single(testCase)
    [u, s, e] = rs3_group_sorted_ids([42]);
    verifyEqual(testCase, u, 42);
    verifyEqual(testCase, s, 1);
    verifyEqual(testCase, e, 1);
end
