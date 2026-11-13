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
#>     solver   package available                                 classes duals
#> 1    highs     highs      TRUE                              LP MILP QP  TRUE
#> 2 clarabel  clarabel      TRUE                        LP QP SOCP CONIC  TRUE
#> 3     scip      scip      TRUE LP MILP QP MIQP QCQP MIQCQP SOCP MISOCP FALSE
#> 4     osqp      osqp      TRUE                                   LP QP  TRUE
#> 5   nloptr    nloptr      TRUE                          NLP QCQP QP LP FALSE
#> 6    ipopt                TRUE                          NLP QCQP QP LP  TRUE
#> 7      roi ROI, slam      TRUE                              LP MILP QP FALSE
#> 8    nlpbb    nloptr      TRUE                       MINLP MIQCQP MIQP FALSE
#>                                                                                    description
#> 1               HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)
#> 2 Clarabel: interior point for LP, convex QP and SOC, PSD, exponential and power cone programs
#> 3         SCIP: branch and bound with spatial branching, global for nonconvex quadratic models
#> 4               OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)
#> 5                            NLopt SLSQP (default): local gradient-based nonlinear programming
#> 6                     IPOPT: interior point for nonlinear programs with exact Hessians (local)
#> 7               ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)
#> 8   NLP branch and bound over NLopt relaxations (global for convex MINLP, heuristic otherwise)
```
