# Logging

The original codebase prints directly with `fprintf`.

This repo adds a tiny logger:

- `rs3_log_init(outdir)` creates `rs3_run.log` in the output folder and sets a default logger.
- `rs3_log(level, fmt, ...)` writes timestamped messages to console + log file.
- `rs3_log_close()` closes the log file.

Example:

```matlab
outdir = fullfile(rs3_repo_root(),'rs3_results','demo');
if ~exist(outdir,'dir'), mkdir(outdir); end
rs3_log_init(outdir,'level','info');
cleanupObj = onCleanup(@() rs3_log_close());

rs3_log('info','Hello from rs3');
rs3_log('debug','This will show only if level is debug');
```

## Suggested refactor pattern

If you want to migrate gradually:

1. Leave inner-loop functions alone for now.
2. Wrap *top-level* stages with `rs3_log('info',...)` and keep `fprintf` inside hot loops.
3. Later, replace `fprintf` in non-hot code paths and keep performance-critical `fprintf` where needed.
