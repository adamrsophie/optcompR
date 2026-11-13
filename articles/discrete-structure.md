# Indicators, SOS, complementarity and logic

``` r

library(optcompR)
```

Some combinatorial structure is easier to state directly than with
hand-written big-M rows: a constraint that holds only when a switch is
on, a set of variables of which only one may be nonzero, or a logical
relation between binary decisions. optcompR has constructors for each
and passes them to solvers either natively or through an exact
reformulation. Nothing is approximated.

| Constructor | Meaning | SCIP | Other MILP solvers |
|----|----|----|----|
| `opt_indicator(m, z, rel, active)` | `rel` holds when `z == active` | native | big-M rows from bounds |
| `opt_sos1(m, x)` | at most one nonzero | native | one binary per element |
| `opt_sos2(m, x)` | at most two adjacent nonzeros | native | one binary per segment |
| `opt_complementarity(m, x, y)` | 0 \le x \perp y \ge 0 | native (SOS1 pairs) | binaries |
| `opt_implies`, `opt_and`, `opt_or` | logic on binaries | linear rows | linear rows |

Models with indicator or SOS constraints are discrete, so they are
classified with the `"MI"` prefix even when no variable is declared
integer.

## Indicator constraints

A machine has a fixed cost of 10 when it runs and a variable cost of 2
per unit; at least 3 units are needed. Production must be zero when the
machine is off, which is an indicator constraint on the complement of
`on`:

``` r

m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 8)
on <- opt_variable(m, "on", type = "binary")
opt_indicator(m, on, x <= 0, active = 0)
opt_constraint(m, x >= 3)
opt_objective(m, 2 * x + 10 * on)
opt_classify(m)$class
#> [1] "MILP"
r <- opt_solve(m, solver = "highs")
opt_value(r, c(x, on))
#> [1] 3 1
r$warnings
#> [1] "Indicator/SOS constraints were reformulated exactly for highs with 0 binary column(s) and 1 row(s)."
```

HiGHS has no indicator constraints, so optcompR rewrote the indicator as
the row x \le 8 \cdot on, computing 8 from the bounds of `x`. The
warning records the reformulation. The big-M value is the tightest one
the bounds allow, which keeps the relaxation as strong as possible.

Without finite bounds there is no valid big-M and the solve stops:

``` r

m2 <- opt_model()
u <- opt_variable(m2, "u", lower = 0)
b <- opt_variable(m2, "b", type = "binary")
opt_indicator(m2, b, u <= 1)
opt_objective(m2, -b)
opt_solve(m2, solver = "highs")
#> Error:
#> ! Indicator constraint 'indicator1' needs finite bounds on its variables to be reformulated exactly for solver 'highs'; add bounds or use solver = "scip".
```

SCIP handles the same model natively:

``` r

opt_objective_value(opt_solve(m2, solver = "scip"))
#> [1] -1
```

## Special ordered sets

An SOS2 set models piecewise linear functions. With breakpoints 0, 1, 2,
3 and weights \lambda that sum to one with at most two adjacent
nonzeros, u = \sum_k k \lambda_k and \sum_k k^2 \lambda_k is the
piecewise linear interpolant of u^2:

``` r

m <- opt_model()
lam <- opt_variable(m, "lambda", dim = 4, lower = 0, upper = 1)
u <- opt_variable(m, "u")
opt_constraint(m, sum(lam) == 1)
opt_constraint(m, u == sum(0:3 * lam))
opt_constraint(m, u >= 1.5)
opt_sos2(m, lam)
opt_objective(m, sum((0:3)^2 * lam))
r <- opt_solve(m)
opt_value(r, lam)
#> [1] 0.0 0.5 0.5 0.0
opt_objective_value(r)  # interpolant at 1.5: (1 + 4) / 2
#> [1] 2.5
```

SOS1 forbids more than one nonzero:

``` r

m <- opt_model("max")
v <- opt_variable(m, "v", dim = 3, lower = 0, upper = c(1, 2, 3))
opt_sos1(m, v)
opt_objective(m, sum(v))
opt_value(opt_solve(m), v)
#> [1] 0 0 3
```

The order of an SOS set is given by `weights` (default `seq_along(x)`);
it matters for SOS2, where only neighbors in that order may be nonzero
together.

## Complementarity

Complementarity constraints 0 \le x_i \perp y_i \ge 0 appear in
equilibrium problems and in the KKT conditions of inner problems.
[`opt_complementarity()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md)
adds one SOS1 pair per element:

``` r

m <- opt_model("max")
a <- opt_variable(m, "a", lower = 0, upper = 5)
b <- opt_variable(m, "b", lower = 0, upper = 5)
opt_constraint(m, a + 2 * b <= 4)
opt_complementarity(m, a, b)
opt_objective(m, a + 1.5 * b)
r <- opt_solve(m)
opt_value(r, c(a, b))
#> [1] 4 0
```

Without the complementarity condition the optimum would mix both.

## Logical constraints

R’s `&&`, `||` and `!` evaluate immediately and cannot build symbolic
constraints. optcompR has explicit functions with exact linear
formulations instead; the negation of a binary `x` is `1 - x`.

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", dim = 3, type = "binary")
all3 <- opt_variable(m, "all3", type = "binary")
any1 <- opt_variable(m, "any1", type = "binary")
opt_and(m, all3, x)   # all3 = x1 and x2 and x3
opt_or(m, any1, x)    # any1 = x1 or x2 or x3
opt_constraint(m, sum(x) == 2)
opt_objective(m, all3 - any1)
r <- opt_solve(m)
opt_value(r, c(all3, any1))
#> [1] 0 1
```

`opt_implies(m, a, b)` adds a \le b. Combined with negation it covers
“at most one of”, “if not a then b” and similar rules:

``` r

m <- opt_model()
open_a <- opt_variable(m, "open_a", type = "binary")
open_b <- opt_variable(m, "open_b", type = "binary")
opt_implies(m, 1 - open_a, open_b)  # if a is closed, b must open
opt_objective(m, 3 * open_a + 2 * open_b)
opt_value(opt_solve(m), c(open_a, open_b))
#> [1] 0 1
```

## Checking solutions

[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
evaluates every indicator and SOS condition at the returned point, in
the original model rather than in the reformulation:

``` r

m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 8)
on <- opt_variable(m, "on", type = "binary")
opt_indicator(m, on, x <= 0, active = 0)
opt_constraint(m, x >= 3)
opt_objective(m, 2 * x + 10 * on)
opt_check(opt_solve(m))$special
#>   constraint      kind violation
#> 1 indicator1 indicator         0
```

## Mixed-integer nonlinear models

The same constructs combine with nonlinear objectives. Mixed-integer
nonlinear models are solved by `"nlpbb"`, a branch and bound over NLopt
relaxations. When the continuous relaxation is convex, pruning by bound
is valid and the search proves global optimality. When it is not known
to be convex, every node is a local solve and the result is labeled
`"locally_optimal"`: a heuristic, not a certificate.

``` r

m <- opt_model()
x <- opt_variable(m, "x", type = "integer", lower = -3, upper = 3)
y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 4)
opt_constraint(m, x + y >= 2.5)
opt_objective(m, exp(0.5 * x) + (y - 1.7)^2)
cls <- opt_classify(m)
c(cls$class, cls$convexity)
#> [1] "MINLP"   "unknown"
r <- opt_solve(m)
r$solver
#> [1] "nlpbb"
opt_status(r)
#> [1] "locally_optimal"
opt_value(r, c(x, y))
#> [1] 1 2
```

This objective is in fact convex, but the convexity analysis does not
prove it for [`exp()`](https://rdrr.io/r/base/Log.html) of an affine
term, so the classification says `"unknown"` and the result is reported
as local. Here it agrees with enumeration of the 35 integer points; in
general a local result needs such a check.
