%RUNTESTS_RS3  Convenience test runner.

rs3_setup
results = runtests('tests');
disp(results)
assertSuccess(results)
