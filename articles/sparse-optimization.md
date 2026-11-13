# Sparse optimization

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Large models are built from vectors and sparse
matrices, not from loops over scalars. Expressions store sparse
triplets, and products with `Matrix` objects stay sparse all the way to
the solver.

## A sparse LP with 20 000 variables

``` r

library(Matrix)
set.seed(1)
n <- 20000
k <- 2000
# every row gets at least one entry so the covering problem is feasible
A <- rsparsematrix(k, n, density = 5 / n, rand.x = function(m) runif(m, 0.5, 2)) +
  sparseMatrix(i = 1:k, j = sample(n, k), x = 1, dims = c(k, n))
b <- runif(k, 1, 10)
cost <- runif(n, 1, 3)

m <- opt_model()
x <- opt_variable(m, "x", dim = n, lower = 0)
build <- system.time({
  opt_constraint(m, A %*% x >= b, name = "cover")
  opt_objective(m, sum(cost * x))
})
build[["elapsed"]]
#> [1] 0.012
sf <- opt_standard_form(m)
dim(sf$A)
#> [1]  2000 20000
class(sf$A)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
r <- opt_solve(m)
opt_status(r)
#> [1] "optimal"
r$timing
#> build solve 
#> 0.001 0.066
```

One vector constraint gives `k` rows; there is one constraint block, so
`opt_dual(r, "cover")` returns all `k` duals.

``` r

summary(opt_dual(r, "cover"))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.2762  0.6998  0.6937  1.0110  2.9710
```

## Sparse quadratic objectives

A least squares problem with a sparse design matrix and a ridge penalty:

``` r

p <- 500
X <- rsparsematrix(2000, p, density = 0.01)
beta <- rnorm(p)
y <- as.vector(X %*% beta) + rnorm(2000, sd = 0.1)
m2 <- opt_model()
bh <- opt_variable(m2, "beta", dim = p)
res <- X %*% bh - y
opt_objective(m2, sum(res^2) + 0.1 * sum(bh^2))
r2 <- opt_solve(m2)
opt_status(r2)
#> [1] "optimal"
sf2 <- opt_standard_form(m2)
c(nonzeros_in_Q = length(sf2$Q@x), size = nrow(sf2$Q))
#> nonzeros_in_Q          size 
#>         23153           500
```

## Guidelines

- Write `A %*% x` and `sum(w * x)` instead of loops when the data is a
  matrix.
- Use one constraint block per family of rows; blocks keep names and
  duals organized and are much faster to build than many scalar
  constraints.
- Parameters avoid rebuilding when only data changes
  ([`vignette("parameterized-models")`](https://adamrsophie.github.io/optcompR/articles/parameterized-models.md)).
