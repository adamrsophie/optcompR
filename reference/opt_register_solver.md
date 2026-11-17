# Register a solver backend

Adds a backend to the registry used by
[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
and
[`opt_solvers()`](https://adamrsophie.github.io/optcompR/reference/opt_solvers.md),
so packages and users can connect other solvers without changing
optcompR. Registering an existing name replaces that backend.

## Usage

``` r
opt_register_solver(
  name,
  solve,
  classes,
  package = character(),
  description = name,
  convex_only = character(),
  global = classes,
  features = list(),
  version = NA_character_,
  priority = 50L
)

opt_unregister_solver(name)
```

## Arguments

- name:

  Backend name used in `opt_solve(solver = )`.

- solve:

  The solve function described below.

- classes:

  Problem classes the backend accepts (see
  [`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)).

- package:

  R packages the backend needs; the backend is available when all of
  them are installed.
  [`character()`](https://rdrr.io/r/base/character.html) for none.

- description:

  One line shown by
  [`opt_solvers()`](https://adamrsophie.github.io/optcompR/reference/opt_solvers.md).

- convex_only:

  Classes accepted only when convexity is proven.

- global:

  Classes the backend solves to global optimality; for other classes it
  is treated as a local solver.

- features:

  Named list of capabilities: `duals`, `ranging`, `warm_start`,
  `time_limit`, `mip_gap`, `iteration_limit` (logical), `cones` (cone
  types such as `"soc"`, `"exp"`, `"pow"`, `"psd"`) and `special`
  (`"indicator"`, `"sos"`: handled natively, otherwise reformulated).

- version:

  A function returning the solver version string, or a string.

- priority:

  Order in automatic selection; lower is tried first. The built-in
  backends use 10 to 90.

## Value

The backend name, invisibly. `opt_unregister_solver()` returns `TRUE`
when a backend was removed.

## The solve function

`solve(sf, cls, control)` receives the
[`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md)
`sf` of the model (parameters substituted; for indicator and SOS
constraints that the backend does not list in `special`, already
reformulated with binary columns and rows), the
[`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)
result `cls`, and `control`, a list with `time_limit` (seconds, `Inf`
for none), `mip_gap` (`NULL` for the solver default), `verbose`,
`warm_start`, `options` (the `options` list given to
[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md))
and `cache`, an environment kept per model and backend in which a
backend may store state, such as a basis, for the next solve of the same
model. `sf$structure_key` changes whenever the model structure does. The
objective is `0.5 x'Qx + c'x + offset` in the sense `sf$sense`, or, when
`sf$objective_kind` is `"nonlinear"`, the lowered expression
`sf$objective_lowered`, which
[`opt_gradient()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md)
and the other derivative functions evaluate. It returns a list with

- `termination`: one of the statuses documented in
  [`opt_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md);

- `primal_status` and `dual_status`: `"feasible_point"`,
  `"infeasible_point"`, `"no_solution"`, `"not_available"` or
  `"unknown"`;

- `x`: the column values (length `length(sf$c)`), `NA` when there is no
  point;

- `objective`: the objective value including `sf$offset`, or `NA`;

- optionally `row_dual` and `col_dual` (shadow prices \\\partial f^\* /
  \partial b\\ of the rows and of the column bounds, in the user's
  sense), `gap`, `bound`, `iterations` (a named numeric vector),
  `message`, `local` (`TRUE` for local solvers, which turns `"optimal"`
  into `"locally_optimal"`), `warm_start`, `warnings` (a character
  vector), `extra` (a list) and `raw` (the solver's own output).

The result is checked: a missing or unknown status, or an `x` of the
wrong length, is an error naming the backend.

## Examples

``` r
# a backend that hands LPs to HiGHS through the highs package
my_lp <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  r <- highs::highs_solve(L = s * sf$c, lower = sf$lower,
                          upper = sf$upper, A = sf$A,
                          lhs = sf$row_lower, rhs = sf$row_upper)
  ok <- r$status == 7L
  list(termination = if (ok) "optimal" else "unknown",
       primal_status = if (ok) "feasible_point" else "no_solution",
       dual_status = "not_available",
       x = if (ok) r$primal_solution else rep(NA_real_, length(sf$c)),
       objective = if (ok) s * r$objective_value + sf$offset else NA_real_)
}
opt_register_solver("my_lp", my_lp, classes = "LP", package = "highs",
                    description = "HiGHS LP through a user backend")
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0, upper = 3)
opt_objective(m, 2 * x)
opt_objective_value(opt_solve(m, solver = "my_lp"))
#> [1] 6
opt_unregister_solver("my_lp")
```
