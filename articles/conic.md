# Conic optimization

``` r

library(optcompR)
```

A conic program minimizes a linear (or convex quadratic) objective over
the intersection of an affine set with convex cones. Many convex
problems that are not linear or quadratic have exact conic forms:
Euclidean norms, geometric means, matrix eigenvalues, entropy and
log-sum-exp. optcompR supports four cones:

| Function | Cone | Constraint |
|----|----|----|
| `opt_soc_constraint(m, x, t)` | second-order | \\x\\\_2 \le t |
| `opt_exp_cone(m, x, y, z)` | exponential | y e^{x/y} \le z,\\ y \> 0 |
| `opt_pow_cone(m, x, y, z, alpha)` | power | x^\alpha y^{1-\alpha} \ge \|z\|,\\ x, y \ge 0 |
| `opt_psd_constraint(m, M)` | semidefinite | M \succeq 0 |

Every argument must be affine in the decision variables once parameters
are substituted. A model with cones is classified as `"SOCP"` when all
its cones are second-order cones and as `"CONIC"` otherwise. Cones are
convex sets, so the model is convex when its objective is linear or a
convex quadratic and its rows are linear; that is recorded in the
classification evidence.

Conic models are solved by [Clarabel](https://clarabel.org), an interior
point solver available from CRAN as the `clarabel` package. SCIP also
accepts second-order cones, which optcompR passes as the equivalent
quadratic rows \\x\\\_2^2 \le t^2,\\ t \ge 0; that path also handles
mixed-integer second-order cone programs (`"MISOCP"`). A solver that
does not support a cone type is refused with a message naming the cone.

## Second-order cones

The distance from the origin to the line x_1 + x_2 = 2 is \sqrt{2}.
Minimize an epigraph variable `t` with \\x\\\_2 \le t:

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
t <- opt_variable(m, "t")
line <- opt_constraint(m, sum(x) == 2, name = "line")
opt_soc_constraint(m, x, t)
opt_objective(m, t)
opt_classify(m)
#> <opt_classification: SOCP, convex>
#>   - syntactic: objective is linear; constraint rows: linear = 1
#>   - structural: 1 cone constraint(s) of type soc; cones are convex sets
#>   - structural: linear objective
r <- opt_solve(m)
r$solver
#> [1] "clarabel"
opt_objective_value(r)
#> [1] 1.414214
opt_value(r, x)
#> [1] 1 1
```

Row duals are shadow prices in the usual optcompR convention, \partial
f^\* / \partial b. Here f^\*(b) = b / \sqrt{2}, so the dual of `line` is
1/\sqrt{2}:

``` r

opt_dual(r, line)
#>      line 
#> 0.7071068
1 / sqrt(2)
#> [1] 0.7071068
```

[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
verifies rows, bounds and cones at the returned point. The cone
residuals appear in `special`. The KKT check is skipped for conic models
because their multipliers live in the dual cones rather than on rows.

``` r

d <- opt_check(r)
d
#> <opt_diagnostics: feasible>
#>   max row violation (scaled):   6.21725e-15
#>   max bound violation (scaled): 0
#>   max integrality violation:    0
#>   max cone/indicator/SOS violation: 4.93957e-09
#>   objective recomputation error: 0
d$special
#>   constraint     kind    violation
#> 1      cone1 cone soc 4.939575e-09
```

### A robust least squares fit

Second-order cones express norms of affine residuals. Minimizing \\Ax -
b\\\_2 exactly:

``` r

set.seed(3)
A <- matrix(rnorm(40), 20)
b <- A %*% c(1, -2) + rnorm(20, sd = 0.1)
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
t <- opt_variable(m, "t")
opt_soc_constraint(m, A %*% x - b, t)
opt_objective(m, t)
r <- opt_solve(m)
opt_value(r, x)
#> [1]  1.018269 -1.993609
qr.solve(A, b)
#>           [,1]
#> [1,]  1.018269
#> [2,] -1.993609
```

The norm of the residual (not its square) is the objective, so the same
solution as ordinary least squares comes out, and the norm itself can
now be bounded, weighted or combined with other norms in one model.

## Semidefinite cones

[`opt_psd_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
takes a square affine matrix expression. Only its upper triangle is
used, so the matrix must be symmetric by construction.
[`opt_sym_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_sym_variable.md)
creates a symmetric matrix of decision variables with k(k+1)/2
underlying variables.

The largest eigenvalue of a symmetric matrix A is the smallest t with
tI - A \succeq 0:

``` r

A <- matrix(c(2, 1, 0,
              1, 3, 1,
              0, 1, 1), 3)
m <- opt_model()
X <- opt_sym_variable(m, "X", 3)
t <- opt_variable(m, "t")
opt_constraint(m, X == t * diag(3) - A)
opt_psd_constraint(m, X)
opt_objective(m, t)
opt_objective_value(opt_solve(m))
#> [1] 3.879385
max(eigen(A)$values)
#> [1] 3.879385
```

A matrix built from expressions works too. The matrix \begin{pmatrix}1 &
x & 0\\ x & 1 & y\\ 0 & y & 1\end{pmatrix} is positive semidefinite
exactly when x^2 + y^2 \le 1, so maximizing x + y gives \sqrt 2:

``` r

m <- opt_model("max")
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
one <- 0 * x + 1
zero <- 0 * x
E <- c(one, x, zero, x, one, y, zero, y, one)
dim(E) <- c(3, 3)
opt_psd_constraint(m, E)
opt_objective(m, x + y)
opt_objective_value(opt_solve(m))
#> [1] 1.414214
```

## Exponential cones

The exponential cone expresses e^x \le t, \log, entropy and log-sum-exp.
Minimizing t subject to e^x \le t and x \ge 1 gives e:

``` r

m <- opt_model()
x <- opt_variable(m, "x", lower = 1)
t <- opt_variable(m, "t")
opt_exp_cone(m, x, 0 * x + 1, t)
opt_objective(m, t)
opt_objective_value(opt_solve(m))
#> [1] 2.718282
exp(1)
#> [1] 2.718282
```

Log-sum-exp \log \sum_i e^{a_i^\top x} \le s is the set of points with
\sum_i u_i \le 1 and e^{a_i^\top x - s} \le u_i, one exponential cone
per term:

``` r

m <- opt_model()
x <- opt_variable(m, "x")
s <- opt_variable(m, "s")
u <- opt_variable(m, "u", dim = 3)
a <- c(1, -1, 0.5)
for (i in 1:3) opt_exp_cone(m, a[i] * x - s, 0 * x + 1, u[i])
opt_constraint(m, sum(u) <= 1)
opt_objective(m, s)
r <- opt_solve(m)
c(x = opt_value(r, x), lse = opt_objective_value(r))
#>          x        lse 
#> -0.2219656  1.0799137
f <- function(x) log(sum(exp(a * x)))
optimize(f, c(-5, 5))
#> $minimum
#> [1] -0.2219163
#> 
#> $objective
#> [1] 1.079914
```

## Power cones

The power cone with exponent \alpha bounds a weighted geometric mean.
With x + y \le 2, the largest z with x^{0.3} y^{0.7} \ge \|z\| is 2
\cdot 0.3^{0.3} 0.7^{0.7}:

``` r

m <- opt_model("max")
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
z <- opt_variable(m, "z")
opt_constraint(m, x + y <= 2)
opt_pow_cone(m, x, y, z, alpha = 0.3)
opt_objective(m, z)
opt_objective_value(opt_solve(m))
#> [1] 1.085763
2 * 0.3^0.3 * 0.7^0.7
#> [1] 1.085763
```

## Mixed-integer second-order cones

With integer variables and only second-order cones the class is
`"MISOCP"`, which SCIP solves to global optimality:

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 2, type = "integer", lower = -5, upper = 5)
t <- opt_variable(m, "t")
opt_soc_constraint(m, x - c(1.4, 2.6), t)
opt_objective(m, t)
opt_classify(m)$class
#> [1] "MISOCP"
r <- opt_solve(m)
opt_value(r, x)
#> [1] 1 3
opt_objective_value(r)
#> [1] 0.5656854
```

## Accuracy

Interior point solvers stop at a relative tolerance near 10^{-8}. The
status `"almost_optimal"` means Clarabel stopped at a reduced accuracy;
check such points with
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
before relying on them. The reported bound is Clarabel’s dual objective,
so `opt_objective_bound(r)` and the objective together bracket the
optimum.
