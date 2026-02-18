function tests = test_rows
%TEST_ROWS  Unit tests for RS3 packed row operations.

tests = functiontests(localfunctions);
end

% ---- rs3_rows_empty ----

function test_empty_default(testCase)
    r = rs3_rows_empty();
    verifyEqual(testCase, double(r.n), 0);
    verifyTrue(testCase, isstruct(r));
    verifyTrue(testCase, isfield(r, 'iSeed'));
    verifyTrue(testCase, isfield(r, 'ix'));
end

function test_empty_preallocated(testCase)
    r = rs3_rows_empty(100);
    verifyEqual(testCase, double(r.n), 0);
    verifyEqual(testCase, numel(r.iSeed), 100);
    verifyEqual(testCase, numel(r.ix), 100);
end

% ---- rs3_rows_count ----

function test_count_struct(testCase)
    r = rs3_rows_empty(50);
    r.n = uint32(25);
    verifyEqual(testCase, rs3_rows_count(r), 25);
end

function test_count_matrix(testCase)
    m = rand(10, 8);
    verifyEqual(testCase, rs3_rows_count(m), 10);
end

function test_count_empty(testCase)
    verifyEqual(testCase, rs3_rows_count(rs3_rows_empty()), 0);
end

function test_count_nonsense(testCase)
    verifyEqual(testCase, rs3_rows_count('hello'), 0);
end

% ---- rs3_rows_vcat ----

function test_vcat_two_empty(testCase)
    a = rs3_rows_empty();
    b = rs3_rows_empty();
    r = rs3_rows_vcat(a, b);
    verifyEqual(testCase, double(r.n), 0);
end

function test_vcat_preserves_data(testCase)
    a = rs3_rows_empty(3);
    a.n = uint32(2);
    a.iSeed(1:2) = uint16([1; 2]);
    a.ix(1:2) = uint16([10; 20]);

    b = rs3_rows_empty(2);
    b.n = uint32(1);
    b.iSeed(1) = uint16(3);
    b.ix(1) = uint16(30);

    r = rs3_rows_vcat(a, b);
    verifyEqual(testCase, double(r.n), 3);
    verifyEqual(testCase, double(r.iSeed(1:3)), [1; 2; 3]);
    verifyEqual(testCase, double(r.ix(1:3)), [10; 20; 30]);
end
