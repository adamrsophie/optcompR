# Continuous problem with integer variables fixed

Fixes every integer and binary variable of a mixed-integer result at its
returned value and solves the remaining continuous problem. Its duals
describe the sensitivity of the continuous part with the discrete
choices held fixed. They are not duals of the mixed-integer problem,
whose optimal value is in general not differentiable in the data; the
returned result is labeled with `fixed_relaxation = TRUE`.

## Usage

``` r
opt_fixed_relaxation(result, ...)
```

## Arguments

- result:

  An `opt_result` of a mixed-integer model with a solution.

- ...:

  Further arguments for
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md).

## Value

An `opt_result` for the continuous problem. The original model is not
changed.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
y <- opt_variable(m, "y", lower = 0)
cap <- opt_constraint(m, x + y <= 4.5)
opt_objective(m, 3 * x + 2 * y)
r <- opt_solve(m)
fr <- opt_fixed_relaxation(r)
opt_dual(fr, cap)
#> c1 
#>  2 
```
