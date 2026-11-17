# Convert a model to a CVXR problem

Builds a CVXR `Problem` from the lowered standard form of a model
(parameters substituted): one CVXR variable vector for the model's
variables, linear and quadratic rows, variable bounds, integrality, and
second order, exponential, power and semidefinite cones. CVXR then
checks the problem against its disciplined convex programming rules, so
a nonconvex quadratic is refused there. Nonlinear expressions other than
quadratics and cones, indicators and SOS constraints have no CVXR
counterpart and give an error.

## Usage

``` r
opt_to_cvxr(model)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

## Value

An `opt_cvxr` list with `problem` (the CVXR `Problem`), `x` (the CVXR
variable, one entry per model variable in the order of `col_names`),
`col_names` and `offset` (the objective constant, which the CVXR
objective includes).

## See also

[`opt_cvxr_values()`](https://adamrsophie.github.io/optcompR/reference/opt_cvxr_values.md)
to read the solution back.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_constraint(m, x[1] + x[2] >= 1)
opt_objective(m, x[1]^2 + x[2]^2)
cv <- opt_to_cvxr(m)
CVXR::psolve(cv$problem)
#> [1] 0.5
opt_cvxr_values(cv)
#> x[1] x[2] 
#>  0.5  0.5 
```
