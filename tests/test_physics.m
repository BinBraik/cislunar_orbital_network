function tests = test_physics
%TEST_PHYSICS  Unit tests for RS3 physics functions.

tests = functiontests(localfunctions);
end

% ---- rs3_core_cr3bp_U_and_derivs ----

function test_potential_at_origin(testCase)
    mu = 0.012150584270572;
    pot = rs3_core_cr3bp_U_and_derivs(0, 0, mu);
    verifyTrue(testCase, isstruct(pot));
    verifyTrue(testCase, isfield(pot, 'U'));
    verifyTrue(testCase, isfinite(pot.U));
end

function test_potential_vectorized(testCase)
    mu = 0.012150584270572;
    x = linspace(-1, 1, 50)';
    y = zeros(50, 1);
    pot = rs3_core_cr3bp_U_and_derivs(x, y, mu);
    verifySize(testCase, pot.U, [50 1]);
end

% ---- rs3_seed_speed ----

function test_seed_speed_empty(testCase)
    v0 = rs3_seed_speed([], 3.13, 0.012);
    verifyEmpty(testCase, v0);
end

function test_seed_speed_positive(testCase)
    mu = 0.012150584270572;
    CJ = 3.13;
    seeds = [0.8 0 pi/4; 1.1 0 -pi/4];
    v0 = rs3_seed_speed(seeds, CJ, mu);
    verifySize(testCase, v0, [2 1]);
    verifyGreaterThanOrEqual(testCase, v0, 0);
end
