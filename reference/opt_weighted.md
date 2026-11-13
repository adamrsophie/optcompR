# Weighted-sum scalarization

Solves the model with the objective \\\sum_k w_k f_k(x)\\ in the model's
sense. With positive weights every solution is Pareto optimal. Weighted
sums cannot reach Pareto points on nonconvex parts of the front (and in
mixed-integer problems they miss unsupported efficient points); use
[`opt_epsilon_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_epsilon_constraint.md)
for those.

## Usage

``` r
opt_weighted(model, objectives, weights, ...)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- objectives:

  A list of scalar expressions.

- weights:

  Nonnegative weights, one per objective.

- ...:

  Further arguments for
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md).

## Value

An `opt_result`, with the individual objective values in
`objective_values`.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 4)
y <- opt_variable(m, "y", lower = 0, upper = 4)
opt_constraint(m, x + y >= 4)
r <- opt_weighted(m, list(cost = 2 * x + y, risk = x + 3 * y), c(1, 1))
r$objective_values
#> cost risk 
#>    8    4 
```
