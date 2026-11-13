# Exact certificates with Z3

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver

Numerical solvers work in floating point with tolerances. For linear
models (LP and MILP) the z3r package gives access to Z3, whose linear
integer and real arithmetic is a decision procedure: its answers `sat`
and `unsat` are exact statements about the model as written.

These functions need z3r and the Z3 library
(`remotes::install_github("adamrsophie/z3r")`).

## Exact feasibility

``` r

m <- opt_model()
x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
opt_constraint(m, 2 * x == 7)
opt_z3_check(m)
#> $status
#> [1] "unsat"
#> 
#> $values
#> NULL
```

No integer solves 2x = 7, and Z3 proves it. A satisfying assignment is a
feasible point, not an optimal one:

``` r

m2 <- opt_model()
x <- opt_variable(m2, "x", type = "integer", lower = 0, upper = 10)
y <- opt_variable(m2, "y", lower = 0)
opt_constraint(m2, 2 * x + 0.1 * y == 7)
opt_objective(m2, y)
opt_z3_check(m2)
#> $status
#> [1] "sat"
#> 
#> $values
#>  x  y 
#>  0 70
opt_z3_check(m2, objective_bound = 9.99)$status
#> [1] "unsat"
```

The second call proves that no feasible point has y \le 9.99: because x
\le 3 is forced, y = 70 - 20x \ge 10.

## Certifying an optimum

[`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)
asks whether any feasible point beats the solver’s objective by more
than a tiny relative margin. `unsat` proves the reported value optimal
for the model.

``` r

k <- opt_model("max")
a <- opt_variable(k, "a", type = "binary")
b <- opt_variable(k, "b", type = "integer", lower = 0, upper = 3)
opt_constraint(k, a + b <= 4)
opt_objective(k, 3 * a + 2 * b)
r <- opt_solve(k)
opt_z3_certify(r)$certified
#> [1] TRUE
```

If the reported value were wrong, Z3 returns a better point as a
counterexample:

``` r

bad <- r
bad$objective_value <- 8
opt_z3_certify(bad)[c("certified", "better_point")]
#> $certified
#> [1] FALSE
#> 
#> $better_point
#> a b 
#> 1 3
```

## Coefficients

Coefficients are doubles. With `coefficients = "decimal"` (the default),
`0.1` means exactly 1/10; with `"binary"` it means the exact binary
value of the double, 3602879701896397 / 2^{55}.

## What Z3 is not used for

Only linear models are translated. Z3’s nonlinear real arithmetic is
incomplete and its semantics differ from numerical optimization, so
quadratic and nonlinear models are refused rather than sent to Z3.
Certification proves the optimal value, not exact feasibility of the
floating point solution vector. Exact solving can be slow for large
integer models.
