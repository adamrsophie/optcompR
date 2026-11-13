# Linear programming

``` r

library(optcompR)
```

A linear program minimizes or maximizes a linear function over a
polyhedron:

\min_x\\ c^\top x \quad \text{subject to}\quad l \le Ax \le u,\\ \\ lb
\le x \le ub.

## Vector variables and matrix data

Expressions are vectorized. A transportation problem with 2 plants and 3
markets uses a 2 by 3 variable and row and column sums.

``` r

cost <- matrix(c(4, 6, 9,
                 5, 3, 7), 2, byrow = TRUE)
supply <- c(35, 50)
demand <- c(20, 30, 25)

m <- opt_model()
ship <- opt_variable(m, "ship", dim = c(2, 3), lower = 0)
for (i in 1:2) opt_constraint(m, sum(ship[i, ]) <= supply[i], name = paste0("supply_", i))
for (j in 1:3) opt_constraint(m, sum(ship[, j]) >= demand[j], name = paste0("demand_", j))
opt_objective(m, sum(cost * ship))
r <- opt_solve(m)
opt_status(r)
#> [1] "optimal"
opt_value(r, ship)
#>      [,1] [,2] [,3]
#> [1,]   20    0    5
#> [2,]    0   30   20
opt_objective_value(r)
#> [1] 355
```

## Ranged and equality constraints

A constraint can have both bounds, either from a relation or from
`lower` and `upper`:

``` r

m2 <- opt_model()
x <- opt_variable(m2, "x", dim = 3, lower = 0)
opt_constraint(m2, sum(x) == 10, name = "total")
opt_constraint(m2, x[1] - x[2], lower = -1, upper = 1, name = "balance")
opt_objective(m2, x[1] + 2 * x[2] + 3 * x[3])
r2 <- opt_solve(m2)
opt_value(r2, x)
#> [1] 5.5 4.5 0.0
opt_constraints(m2)
#>      name rows   kind
#> 1   total    1 affine
#> 2 balance    1 affine
```

## The standard form

[`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md)
shows exactly what the solver receives: constants are moved to the
bounds and every row is l \le a^\top x \le u.

``` r

sf <- opt_standard_form(m2)
sf$A
#> 2 x 3 sparse Matrix of class "dgCMatrix"
#>            
#> [1,] 1  1 1
#> [2,] 1 -1 .
cbind(lower = sf$row_lower, upper = sf$row_upper)
#>      lower upper
#> [1,]    10    10
#> [2,]    -1     1
```

## Reduced costs and sensitivity ranges

``` r

opt_reduced_cost(r2, x)
#> x[1] x[2] x[3] 
#>  0.0  0.0  1.5
s <- opt_sensitivity(r2)
s$costs
#>   variable coefficient lower upper objective_at_lower objective_at_upper
#> 1     x[1]           1  -Inf     2               -Inf               20.0
#> 2     x[2]           2   1.0     5               10.0               28.0
#> 3     x[3]           3   1.5   Inf               14.5               14.5
s$rows
#>       row row_lower row_upper dual bound_lower bound_upper objective_at_lower
#> 1   total        10        10  1.5           1         Inf                1.0
#> 2 balance        -1         1 -0.5          -1          10               15.5
#>   objective_at_upper
#> 1                Inf
#> 2                 10
```

The cost ranges say how far each objective coefficient can move before
the optimal basis changes; the row ranges say the same for the bounds.

## Export

LP and MILP models can be written to MPS or LP files for other tools:

``` r

f <- tempfile(fileext = ".mps")
opt_write_model(m2, f)
#> WARNING: Column names are not present, or contain duplicates: using names with prefix "c", beginning with suffix 0
#> WARNING: Row    names are not present, or contain duplicates: using names with prefix "r", beginning with suffix 0
#> Writing the model to /tmp/Rtmp1hgZu4/file1eca63e82e73.mps
head(readLines(f))
#> [1] "NAME        " "ROWS"         " N  Obj     " " E  r0      " " L  r1      "
#> [6] "COLUMNS"
```
