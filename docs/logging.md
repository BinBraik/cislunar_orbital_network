# Logging

The codebase prints progress directly with `fprintf` using the prefixes `[atlas]`,
`[overlap]`, and `[traj]` to indicate which pipeline stage is active.

Verbose output is controlled by `cfg.io.verbose` (default: `true`).
Progress reporting inside PARFOR loops is controlled by `cfg.diag.progress` and
`cfg.par.progress_every`.

To suppress all output, set both to `false` before running a runner script:

```matlab
cfg = atlas_cfg_defaults();
cfg.io.verbose      = false;
cfg.diag.progress   = false;
```

## Adding a custom logger

If you want to redirect output to a log file, wrap the runner in a MATLAB diary
or replace `fprintf` calls with your own logging backend. The pattern used in
runner scripts is straightforward to intercept.
