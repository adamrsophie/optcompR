# Exact feasibility check with Z3

Translates the constraints and bounds of a linear model (LP or MILP,
after parameter substitution) to linear integer and real arithmetic and
asks Z3 whether a feasible point exists. For these theories Z3 is a
decision procedure, so `"sat"` and `"unsat"` are exact answers about the
model as written, without floating-point tolerances. A satisfying
assignment is a feasible point, not an optimal one.

## Usage

``` r
opt_z3_check(
  model,
  objective_bound = NULL,
  coefficients = c("decimal", "binary")
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md)
  with linear constraints.

- objective_bound:

  Optional number: also require the objective to be at most
  (minimization) or at least (maximization) this value.

- coefficients:

  `"decimal"` or `"binary"`, see Details.

## Value

A list with `status` (`"sat"`, `"unsat"` or `"unknown"`) and, when
satisfiable, `values`, a named feasible point.

## Details

Coefficients and bounds are doubles. With `coefficients = "decimal"`
each is read as its 15-significant-digit decimal (so `0.1` means 1/10);
with `"binary"` as the exact binary value of the double.

## See also

[`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
opt_constraint(m, 2 * x == 7)
opt_z3_check(m)$status
#> [1] "unsat"
```
