# Contributing to optcompR

- Open an issue describing the bug or feature before a large change.
- Work on a branch and open a pull request against `main`.
- Every new feature needs a working computation and a test that checks
  it against an independent reference (an analytic value, another solver
  or a finite difference).
- Run `devtools::document()` and keep `R CMD check` free of warnings; CI
  runs it on macOS, Windows and Ubuntu.
- New solver backends register through `register_solver()` with their
  classes, convexity requirements and features.
