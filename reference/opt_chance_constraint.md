# Chance constraints

Requires \\P(a^\top x \le b) \ge p\\ for a random coefficient vector
\\a\\.

## Usage

``` r
opt_chance_constraint(
  model,
  x,
  b,
  prob,
  mean = NULL,
  Sigma = NULL,
  samples = NULL,
  method = c("gaussian", "scenario"),
  name = NULL
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  Affine vector expression, usually decision variables.

- b:

  Right-hand side.

- prob:

  Required probability.

- mean, Sigma:

  Mean vector and covariance matrix (Gaussian method).

- samples:

  Matrix of sampled coefficient vectors, one per row (scenario method).

- method:

  `"gaussian"` or `"scenario"`.

- name:

  Base name of the generated constraints.

## Value

The constraint reference, invisibly.

## Details

- `method = "gaussian"`: \\a \sim N(\mu, \Sigma)\\. The constraint is
  equivalent to \\\mu^\top x + \Phi^{-1}(p) \\\Sigma^{1/2} x\\\_2 \le
  b\\, a second-order cone for \\p \ge 0.5\\. This is exact.

- `method = "scenario"`: the constraint is imposed for every row of
  `samples`. This is the scenario approximation: with \\N\\ independent
  samples and \\d\\ decision variables, the solution violates the chance
  constraint with probability at most \\\epsilon\\ when \\N \ge
  \frac{2}{\epsilon}(\ln\frac{1}{\beta} + d)\\ (Calafiore and Campi),
  with confidence \\1 - \beta\\. It is labeled an approximation, not an
  exact reformulation.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_chance_constraint(m, x, b = 10, prob = 0.95, mean = c(1, 1),
                      Sigma = diag(c(0.04, 0.01)))
opt_objective(m, sum(x))
if (requireNamespace("clarabel", quietly = TRUE)) opt_value(opt_solve(m), x)
#> [1] 1.743499 6.973984
```
