# Conic constraints

Convex cone constraints on affine expressions, solved by conic solvers
(Clarabel) and, for second-order cones, by SCIP:

## Usage

``` r
opt_soc_constraint(model, x, t, name = NULL)

opt_exp_cone(model, x, y, z, name = NULL)

opt_pow_cone(model, x, y, z, alpha, name = NULL)

opt_psd_constraint(model, M, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x, y, z, t:

  Affine expressions; `x` in `opt_soc_constraint()` may be a vector, all
  others are scalars.

- name:

  Optional name.

- alpha:

  Exponent of the power cone.

- M:

  Square affine matrix expression.

## Value

An `opt_special_constraint` reference, invisibly.

## Details

- `opt_soc_constraint()`: \\\\x\\\_2 \le t\\ (second-order cone).

- `opt_exp_cone()`: \\y \exp(x / y) \le z\\ with \\y \> 0\\ (the closure
  of the exponential cone, as in Clarabel and SCS).

- `opt_pow_cone()`: \\x^\alpha y^{1-\alpha} \ge \|z\|\\ with \\x, y \ge
  0\\ and \\0 \< \alpha \< 1\\.

- `opt_psd_constraint()`: the symmetric matrix `M` is positive
  semidefinite. Only the upper triangle of `M` is used, so symmetry must
  hold by construction (for example with
  [`opt_sym_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_sym_variable.md)).

All arguments must be affine in the decision variables after parameter
substitution. Models with cones are classified as `"SOCP"` (second-order
cones only) or `"CONIC"`, with `"MI"` prefixed when integer variables
are present; they are convex when the objective is linear or a convex
quadratic and all rows are linear.

## Examples

``` r
# distance from the origin to the line x + y = 2
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
t <- opt_variable(m, "t")
opt_constraint(m, sum(x) == 2)
opt_soc_constraint(m, x, t)
opt_objective(m, t)
opt_objective_value(opt_solve(m))  # sqrt(2)
#> [1] 1.414214
```
