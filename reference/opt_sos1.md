# Special ordered sets and complementarity

`opt_sos1()`: at most one element of `x` is nonzero. `opt_sos2()`: at
most two elements are nonzero and they are adjacent in the order given
by `weights`. `opt_complementarity()`: \\x_i y_i = 0\\ elementwise for
nonnegative variables, the complementarity condition \\0 \le x \perp y
\ge 0\\, as SOS1 pairs.

## Usage

``` r
opt_sos1(model, x, weights = NULL, name = NULL)

opt_sos2(model, x, weights = NULL, name = NULL)

opt_complementarity(model, x, y, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x, y:

  Decision variables (a block or a slice of one).

- weights:

  Distinct weights that order the set; defaults to `seq_along(x)`.

- name:

  Optional name.

## Value

An `opt_special_constraint` reference (a list of them for
`opt_complementarity()`), invisibly.

## Details

SCIP handles these natively. For other mixed-integer solvers they are
reformulated exactly with one binary per element (SOS1) or segment
(SOS2), which requires finite bounds.

## Examples

``` r
# piecewise linear interpolation of f(u) = u^2 on the grid 0, 1, 2, 3
m <- opt_model()
lam <- opt_variable(m, "lambda", dim = 4, lower = 0, upper = 1)
u <- opt_variable(m, "u")
opt_constraint(m, sum(lam) == 1)
opt_constraint(m, u == sum(0:3 * lam))
opt_constraint(m, u >= 1.5)
opt_sos2(m, lam)
opt_objective(m, sum((0:3)^2 * lam))
opt_objective_value(opt_solve(m))  # 2.5, the interpolant at 1.5
#> [1] 2.5
```
