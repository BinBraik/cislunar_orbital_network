%RUNTESTS  Convenience test runner.

setup
results = runtests('tests');
disp(results)
assertSuccess(results)
