# Linear programming sensitivity report

For an optimal linear program solved by HiGHS, reports the ranges over
which objective coefficients and row bounds can move while the optimal
basis stays optimal, together with the objective value at the ends of
each range. Ranges are local to the returned basis; at degenerate optima
different bases give different ranges, and outside a range the optimal
value is piecewise linear rather than following the dual.

## Usage

``` r
opt_sensitivity(result)
```

## Arguments

- result:

  An `opt_result` for an LP solved with `solver = "highs"`.

## Value

A list of two data frames: `costs` (variable, coefficient, lowest and
highest coefficient keeping the basis optimal, and the objective at
those values) and `rows` (row, active bound, dual, lowest and highest
bound value and the objective there).

## Details

Row bounds are reported in standard-form terms
(`row_lower <= a'x <= row_upper`, constants moved to the bounds), so for
`a'x + c <= b` the bound is `b - c`.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
opt_constraint(m, x + y <= 4)
opt_constraint(m, x <= 2)
opt_objective(m, 3 * x + 2 * y)
opt_sensitivity(opt_solve(m))
#> $costs
#>   variable coefficient lower upper objective_at_lower objective_at_upper
#> 1        x           3     2   Inf                  8                Inf
#> 2        y           2     0     3                  6                 12
#> 
#> $rows
#>   row row_lower row_upper dual bound_lower bound_upper objective_at_lower
#> 1  c1      -Inf         4    2           2         Inf                  6
#> 2  c2      -Inf         2    1           0           4                  8
#>   objective_at_upper
#> 1                Inf
#> 2                 12
#> 
```
