# Constraint values, slacks and activity

Evaluates each row of a constraint block from the original expression at
the solution, independently of the standard form.

## Usage

``` r
opt_slack(result, constraint, tol = 1e-07)
```

## Arguments

- result:

  An `opt_result`.

- constraint:

  An `opt_constraint` reference or block name.

- tol:

  Absolute tolerance, scaled by the bound magnitude, for marking rows
  active.

## Value

A data frame with the row value, lower and upper bound, slack (distance
to the nearest bound, negative when violated) and an `active` flag.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
cn <- opt_constraint(m, x <= 3)
opt_objective(m, x)
opt_slack(opt_solve(m), cn)
#>   row value lower upper slack active
#> 1  c1     0  -Inf     0     0   TRUE
```
