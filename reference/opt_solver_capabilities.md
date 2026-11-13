# Capabilities of one solver backend

Capabilities of one solver backend

## Usage

``` r
opt_solver_capabilities(solver)
```

## Arguments

- solver:

  Backend name, see
  [`opt_solvers()`](https://adamrsophie.github.io/optcompR/reference/opt_solvers.md).

## Value

A list with the problem classes, the classes accepted only when convex,
the classes solved to global optimality, feature flags (duals, ranging,
warm starts, time limits, MIP gaps, iteration limits), availability and
the installed version.

## Examples

``` r
opt_solver_capabilities("highs")
#> $name
#> [1] "highs"
#> 
#> $package
#> [1] "highs"
#> 
#> $available
#> [1] TRUE
#> 
#> $version
#> [1] "1.15.1.2"
#> 
#> $classes
#> [1] "LP"   "MILP" "QP"  
#> 
#> $convex_only
#> [1] "QP"
#> 
#> $global
#> [1] "LP"   "MILP" "QP"  
#> 
#> $features
#> $features$duals
#> [1] TRUE
#> 
#> $features$ranging
#> [1] TRUE
#> 
#> $features$warm_start
#> [1] TRUE
#> 
#> $features$time_limit
#> [1] TRUE
#> 
#> $features$mip_gap
#> [1] TRUE
#> 
#> $features$iteration_limit
#> [1] TRUE
#> 
#> 
#> $description
#> [1] "HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)"
#> 
```
