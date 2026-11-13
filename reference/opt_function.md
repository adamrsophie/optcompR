# Wrap an R function as a nonlinear expression

Creates a scalar expression \\f(x)\\ from an R function of a numeric
vector, for objectives or constraints that are not written with optcompR
operators: simulation-based losses, functions written for jaxR, or
functions of a model from modelcompR. Solvers need a gradient, and its
source is explicit:

## Usage

``` r
opt_function(fn, x, gradient, hessian = NULL, name = "f", step = 1e-06)
```

## Arguments

- fn:

  Function of a numeric vector returning a number.

- x:

  Expression the function is applied to (for example a variable vector).

- gradient:

  A gradient function, `"jaxR"` or `"finite_difference"`.

- hessian:

  Optional function returning the Hessian matrix.

- name:

  Label used in printing.

- step:

  Relative step for finite differences.

## Value

A scalar nonlinear `opt_expr`.

## Details

- a function returning the gradient vector (exact, as supplied);

- `"jaxR"`: reverse-mode automatic differentiation with
  [`jaxR::grad()`](https://adamrsophie.github.io/jaxR/reference/grad.html).
  `fn` must be traceable by jaxR and receives a jaxR array, whose
  indexing is 0-based as in JAX (`v[0]` is the first element), so the
  same `fn` is evaluated by jaxR for values too. It runs with 64-bit
  types enabled (`jax_enable_x64`);

- `"finite_difference"`: central differences with relative step `step`.
  This is an approximation with error of order `step^2` (plus rounding
  of order `eps / step`); results record it, and it is never chosen
  implicitly.

Hessians are available only when `hessian` is a function.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, start = c(-1, 1))
rosen <- function(v) (1 - v[1])^2 + 100 * (v[2] - v[1]^2)^2
rosen_grad <- function(v) c(-2 * (1 - v[1]) - 400 * v[1] * (v[2] - v[1]^2),
                            200 * (v[2] - v[1]^2))
opt_objective(m, opt_function(rosen, x, gradient = rosen_grad,
                              name = "rosenbrock"))
if (requireNamespace("nloptr", quietly = TRUE)) opt_value(opt_solve(m), x)
#> [1] 1 1
```
