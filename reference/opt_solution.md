# Solution table

Solution table

## Usage

``` r
opt_solution(result)
```

## Arguments

- result:

  An `opt_result`.

## Value

A data frame with one row per variable: name, value, bounds, type and
(for continuous problems) reduced cost.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = c(1, 2))
opt_objective(m, sum(x))
opt_solution(opt_solve(m))
#>   name value lower upper       type reduced_cost
#> 1 x[1]     1     1   Inf continuous            1
#> 2 x[2]     2     2   Inf continuous            1
```
