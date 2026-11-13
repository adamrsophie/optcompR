# Objective value of a result

The objective value is recomputed from the model's objective expression
at the returned point, including constant terms, so it does not depend
on how a backend handles offsets. The value reported by the solver is
kept as `solver = TRUE`;
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
compares the two. For mixed-integer problems `opt_objective_bound()`
returns the best proven bound and `opt_gap()` the relative gap.

## Usage

``` r
opt_objective_value(result, solver = FALSE)

opt_objective_bound(result)

opt_gap(result)
```

## Arguments

- result:

  An `opt_result`.

- solver:

  Return the solver-reported value instead.

## Value

A number (`NA` when there is no primal point).

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 2)
opt_objective(m, 3 * x + 10)
opt_objective_value(opt_solve(m))
#> [1] 16
```
