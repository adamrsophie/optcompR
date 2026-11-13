# Certify a MILP or LP optimum with Z3

Checks with exact arithmetic that no feasible point has an objective
better than the reported optimum by more than `tol * max(1, |f*|)`. An
`"unsat"` answer is a proof (for the model as translated, see
[`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md))
that the solver's objective value is optimal up to that margin. It does
not prove that the returned point itself is feasible in exact
arithmetic; floating-point solutions satisfy constraints within solver
tolerances. Exact solving can be slow for large integer models.

## Usage

``` r
opt_z3_certify(result, tol = 1e-09, coefficients = c("decimal", "binary"))
```

## Arguments

- result:

  An `opt_result` of an LP or MILP with a solution.

- tol:

  Relative margin.

- coefficients:

  See
  [`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md).

## Value

A list with `certified` (`TRUE` when no feasible point beats the
reported value by more than the margin and some feasible point attains
it within the margin, `FALSE` otherwise, `NA` when Z3 returns unknown),
`attained`, `status`, the `threshold` that was tested and, when a better
point exists, `better_point`.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", type = "binary")
y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 3)
opt_constraint(m, x + y <= 4)
opt_objective(m, 3 * x + 2 * y)
opt_z3_certify(opt_solve(m))$certified
#> [1] TRUE
```
