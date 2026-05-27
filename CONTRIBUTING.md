# Contributing

Thanks for contributing to **cr3bp_reachable_sets**.

## Development workflow

1. Create a feature branch.
2. Keep changes scoped and focused.
3. Run tests before opening a PR.
4. Update docs when behavior or interfaces change.

## MATLAB style guidelines

- Prefer small, composable functions in `src/`.
- Use descriptive names with `rs3_` / `rs4_` prefixes for pipeline functions.
- Validate key inputs at boundaries (`*_validate` style) when adding new entry points.
- Keep runner scripts in `scripts/`; avoid hard-coding machine-specific paths.

## Testing

Use the provided test runner:

```matlab
run('tests/runtests.m')
```

Add/update tests in `tests/` for any behavior changes.

## Documentation

If your change affects architecture or logging expectations, update:

- `docs/architecture.md`
- `docs/logging.md`

Also update `README.md` when setup or execution instructions change.
