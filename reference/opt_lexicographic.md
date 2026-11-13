# Lexicographic optimization

Optimizes the objectives in priority order. After objective \\k\\ is
optimized with value \\f_k^\*\\, the constraint \\f_k(x) \le f_k^\* +
\mathrm{tol} \max(1, \|f_k^\*\|)\\ (or \\\ge\\ for maximization) is
added before the next objective. The relative tolerance keeps later
stages feasible despite solver tolerances; it also means later stages
may trade a tiny amount of earlier objectives.

## Usage

``` r
opt_lexicographic(model, objectives, tol = 1e-07, ...)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- objectives:

  A list of scalar expressions.

- tol:

  Relative tolerance on earlier objectives.

- ...:

  Further arguments for
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md).

## Value

A list of `opt_result` objects, one per stage; the last one is the
lexicographic solution. Objective values of all stages are in its
`objective_values`.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 4)
y <- opt_variable(m, "y", lower = 0, upper = 4)
opt_constraint(m, x + y >= 4)
st <- opt_lexicographic(m, list(cost = x + y, risk = x))
opt_value(st[[2]], c(x, y))
#> [1] 0 4
```
