# Robust linear constraints

Imposes \\a^\top x \le b\\ for every coefficient vector \\a\\ in an
uncertainty set, through an exact reformulation:

## Usage

``` r
opt_robust_constraint(
  model,
  x,
  a0,
  b,
  uncertainty = c("box", "ellipsoid", "budget"),
  delta = NULL,
  P = NULL,
  gamma = NULL,
  name = NULL
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  Affine vector expression, usually decision variables.

- a0:

  Nominal coefficients.

- b:

  Right-hand side.

- uncertainty:

  `"box"`, `"ellipsoid"` or `"budget"`.

- delta:

  Half-widths of the box (box and budget sets).

- P:

  Matrix with `length(x)` rows shaping the ellipsoid.

- gamma:

  Budget of uncertainty, between 0 and `length(x)`.

- name:

  Base name of the generated constraints.

## Value

A list of the generated constraint references, invisibly.

## Details

- `"box"`: \\a_i \in \[a\_{0i} - \delta_i, a\_{0i} + \delta_i\]\\,
  giving \\a_0^\top x + \delta^\top \|x\| \le b\\ (linear, with
  auxiliary variables for \\\|x\|\\).

- `"ellipsoid"`: \\a = a_0 + P u\\ with \\\\u\\\_2 \le 1\\, giving
  \\a_0^\top x + \\P^\top x\\\_2 \le b\\ (a second-order cone).

- `"budget"`: box deviations of which at most `gamma` act at once
  (Bertsimas and Sim), reformulated by linear duality.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, delta = c(0.5, 0.2))
opt_objective(m, 3 * x[1] + 2 * x[2])
opt_value(opt_solve(m), x)
#> [1] 6.666667 0.000000
```
