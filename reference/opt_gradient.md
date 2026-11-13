# Exact derivatives of the objective and constraints

`opt_gradient()` returns the gradient of the objective, `opt_jacobian()`
the sparse Jacobian of all constraint rows (row order of
[`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md)),
and `opt_hessian()` the Hessian of the objective, or, when `multipliers`
are given, the Hessian of the Lagrangian \\L = f(x) - \sum_i y_i
g_i(x)\\ in the sign convention of
[`opt_dual()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md).
Derivatives are exact: polynomial parts are differentiated structurally
(for \\f = \frac12 x^\top Q x + c^\top x\\ the gradient is \\\frac12
(Q + Q^\top) x + c\\), and nonlinear nodes use the chain rule with
sparse forward-mode Jacobians and exact second derivatives. No finite
differences are used. Parameters take their current values. Functions
that are not differentiable at a point (for example
[`abs()`](https://rdrr.io/r/base/MathFun.html) at zero) return a
one-sided value; check the point when it matters.

## Usage

``` r
opt_gradient(model, at)

opt_jacobian(model, at)

opt_hessian(model, at, multipliers = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- at:

  Point: a numeric vector with one value per variable (named by variable
  names, or in model order), or an `opt_result`, in which case the
  parameter values stored in the result are used too.

- multipliers:

  Optional row multipliers for the Lagrangian Hessian.

## Value

A named numeric vector (gradient) or a sparse `Matrix` (Jacobian,
Hessian) with dimnames.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
opt_objective(m, x^2 + 3 * x * y + exp(y))
opt_constraint(m, x * y <= 4)
opt_gradient(m, c(1, 2))
#>        x        y 
#>  8.00000 10.38906 
opt_jacobian(m, c(1, 2))
#> 1 x 2 sparse Matrix of class "dgCMatrix"
#>    x y
#> c1 2 1
opt_hessian(m, c(1, 2))
#> 2 x 2 sparse Matrix of class "dsCMatrix"
#>   x        y
#> x 2 3.000000
#> y 3 7.389056
```
