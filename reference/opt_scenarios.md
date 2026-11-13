# Solve a model over parameter scenarios

Sets a parameter to each value in turn, solves, and collects the status
and objective. The parameter is restored afterwards. Each solve rebuilds
the numerical data from the cached structure.

## Usage

``` r
opt_scenarios(model, parameter, values, ..., keep = FALSE)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- parameter:

  A parameter expression or name.

- values:

  A list of values or a numeric vector (one scenario per element for
  scalar parameters).

- ...:

  Further arguments for
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md).

- keep:

  Keep the full results in the `results` attribute.

## Value

A data frame with one row per scenario: index, status and objective
value.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
d <- opt_parameter(m, "d", value = 1)
opt_constraint(m, x >= d)
opt_objective(m, 2 * x)
opt_scenarios(m, d, c(1, 2, 3))
#>   scenario value  status objective
#> 1        1     1 optimal         2
#> 2        2     2 optimal         4
#> 3        3     3 optimal         6
```
