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
#>     solver   package available
#> 1    highs     highs      TRUE
#> 2 clarabel  clarabel      TRUE
#> 3     scip      scip      TRUE
#> 4     osqp      osqp      TRUE
#> 5      scs       scs      TRUE
#> 6   nloptr    nloptr      TRUE
#> 7    ipopt                TRUE
#> 8      roi ROI, slam      TRUE
#> 9    nlpbb    nloptr      TRUE
#>                                             classes duals
#> 1                                        LP MILP QP  TRUE
#> 2                                  LP QP SOCP CONIC  TRUE
#> 3 LP MILP QP MIQP QCQP MIQCQP SOCP MISOCP NLP MINLP FALSE
#> 4                                             LP QP  TRUE
#> 5                                  LP QP SOCP CONIC  TRUE
#> 6                                    NLP QCQP QP LP FALSE
#> 7                                    NLP QCQP QP LP  TRUE
#> 8                                        LP MILP QP FALSE
#> 9                                 MINLP MIQCQP MIQP FALSE
#>                                                                                                                    description
#> 1                                               HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)
#> 2                                 Clarabel: interior point for LP, convex QP and SOC, PSD, exponential and power cone programs
#> 3                          SCIP: branch and bound with spatial branching, global for nonconvex quadratic and polynomial models
#> 4                                               OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)
#> 5 SCS: operator splitting for LP, convex QP and SOC, PSD, exponential and power cone programs (first-order, moderate accuracy)
#> 6                                                            NLopt SLSQP (default): local gradient-based nonlinear programming
#> 7                                                     IPOPT: interior point for nonlinear programs with exact Hessians (local)
#> 8                                               ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)
#> 9                                   NLP branch and bound over NLopt relaxations (global for convex MINLP, heuristic otherwise)
```
