# Derivatives of the optimal solution with respect to parameters

Differentiates the optimal point \\x^\*(p)\\ and optimal value
\\f^\*(p)\\ of a continuous model with respect to parameters, by the
implicit function theorem applied to the KKT conditions on the active
set at the solution. With \\H\\ the Hessian of the Lagrangian \\L = f -
y^\top g\\ in the free variables (those not at a bound), \\J_A\\ the
Jacobian of the active rows and \\y_A\\ their multipliers, the
derivatives solve \$\$\begin{pmatrix} H & -J_A^\top \\ J_A & 0
\end{pmatrix} \begin{pmatrix} dx/dp \\ dy_A/dp \end{pmatrix} =
-\begin{pmatrix} \nabla^2\_{xp} L \\ \partial g_A / \partial p
\end{pmatrix}\$\$ and \\df^\*/dp = \partial L / \partial p\\ (envelope
theorem). All second derivatives are exact.

## Usage

``` r
opt_solution_jacobian(result, parameter, tol = 1e-07)
```

## Arguments

- result:

  A solved `opt_result` with a primal solution.

- parameter:

  A parameter expression or name.

- tol:

  Tolerance for active rows, bounds and zero multipliers.

## Value

A list with `solution` (matrix: variables by parameter elements),
`objective` (gradient of the optimal value), `active` (names of active
rows and bounds) and `multipliers` (`"solver"` or `"estimated"`).

## Details

The result is valid when the solution is a strict local optimum that
satisfies linear independence of the active constraint gradients (LICQ),
the second-order sufficient condition and strict complementarity. The
function stops when the KKT matrix is singular (LICQ or second-order
sufficiency fails) and warns when an active row or bound has a zero
multiplier (strict complementarity fails, so the solution map may only
have one-sided derivatives). The active set is held fixed, so the
derivatives are local: a large parameter change can change the active
set. Mixed-integer and conic models are refused; for a mixed-integer
model, differentiate the
[`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
instead.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
p <- opt_parameter(m, "p", value = 1)
opt_constraint(m, sum(x) <= p, name = "budget")
opt_objective(m, sum((x - c(2, 1))^2))
r <- opt_solve(m)
opt_solution_jacobian(r, p)  # dx/dp = (0.5, 0.5)
#> $solution
#>        p
#> x[1] 0.5
#> x[2] 0.5
#> 
#> $objective
#>  p 
#> -2 
#> 
#> $active
#> [1] "budget"
#> 
#> $multipliers
#> [1] "solver"
#> 
```
