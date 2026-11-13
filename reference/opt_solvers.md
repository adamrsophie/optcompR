# Available optimization solvers

Lists the registered backends, whether their R packages are installed,
and the problem classes each one accepts. Capabilities describe what the
R binding exposes, which can be less than the native solver supports.

## Usage

``` r
opt_solvers()
```

## Value

A data frame with one row per backend.

## Examples

``` r
opt_solvers()
#>   solver   package available                     classes duals
#> 1  highs     highs      TRUE                  LP MILP QP  TRUE
#> 2   scip      scip      TRUE LP MILP QP MIQP QCQP MIQCQP FALSE
#> 3   osqp      osqp      TRUE                       LP QP  TRUE
#> 4 nloptr    nloptr      TRUE              NLP QCQP QP LP FALSE
#> 5    roi ROI, slam      TRUE                  LP MILP QP FALSE
#>                                                                            description
#> 1       HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)
#> 2 SCIP: branch and bound with spatial branching, global for nonconvex quadratic models
#> 3       OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)
#> 4                    NLopt SLSQP (default): local gradient-based nonlinear programming
#> 5       ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)
```
