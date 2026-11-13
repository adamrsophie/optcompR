# Dual values and reduced costs

Duals follow one convention for every backend and both senses: the dual
of a constraint row is the shadow price \\\partial f^\* / \partial b\\,
the rate of change of the optimal objective when the active bound of the
row is increased. The reduced cost of a variable is the same derivative
for its active bound. Under this convention the Lagrangian is \\L(x, y,
z) = f(x) - \sum_i y_i (g_i(x) - b_i) - \sum_j z_j (x_j - l_j)\\ and
stationarity reads \\\nabla f(x) - J(x)^\top y - z = 0\\. For a
minimization, rows bounded above have \\y \le 0\\ and rows bounded below
\\y \ge 0\\; signs flip for maximization.

## Usage

``` r
opt_dual(result, constraint)

opt_reduced_cost(result, x)
```

## Arguments

- result:

  An `opt_result`.

- constraint:

  An `opt_constraint` reference or block name.

- x:

  A variable expression or block name.

## Value

A numeric vector named by row.

## Details

Duals are defined for continuous problems. For mixed-integer problems
they are not available, because the optimal value is not differentiable
in the data in general;
[`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
builds the continuous problem with integers fixed, whose duals are a
labeled local quantity.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
cap <- opt_constraint(m, x + y <= 4)
opt_constraint(m, x <= 2)
opt_objective(m, 3 * x + 2 * y)
r <- opt_solve(m)
opt_dual(r, cap)
#> c1 
#>  2 
opt_reduced_cost(r, x)
#> x 
#> 0 
```
