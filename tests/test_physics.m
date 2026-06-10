function tests = test_physics
%TEST_PHYSICS  Unit tests for RS3 physics functions.

tests = functiontests(localfunctions);
end

% ---- cr3bp_potential ----

function test_potential_at_origin(testCase)
    mu = 0.012150584270572;
    pot = cr3bp_potential(0, 0, mu);
    verifyTrue(testCase, isstruct(pot));
    verifyTrue(testCase, isfield(pot, 'U'));
    verifyTrue(testCase, isfinite(pot.U));
end

function test_potential_vectorized(testCase)
    mu = 0.012150584270572;
    x = linspace(-1, 1, 50)';
    y = zeros(50, 1);
    pot = cr3bp_potential(x, y, mu);
    verifySize(testCase, pot.U, [50 1]);
end

