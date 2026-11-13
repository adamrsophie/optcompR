# Epsilon-constraint method for two objectives

Traces a Pareto front by optimizing `primary` subject to
`secondary <= epsilon` (for minimization; `>=` for maximization) over a
grid of `epsilon` values. Unlike weighted sums, this reaches points on
nonconvex parts of the front. Points can be weakly Pareto optimal; a
second lexicographic stage would remove dominated ties.

## Usage

``` r
opt_epsilon_constraint(model, primary, secondary, epsilon, ...)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- primary, secondary:

  Scalar expressions.

- epsilon:

  Numeric vector of bounds on `secondary`.

- ...:

  Further arguments for
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md).

## Value

A data frame with `epsilon`, status, the two objective values, with the
results in the `results` attribute.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 4)
y <- opt_variable(m, "y", lower = 0, upper = 4)
opt_constraint(m, x + y >= 4)
opt_epsilon_constraint(m, 2 * x + y, x + 3 * y, epsilon = c(4, 8, 12))
#>   epsilon  status primary secondary
#> 1       4 optimal       8         4
#> 2       8 optimal       6         8
#> 3      12 optimal       4        12
```
