# Benchmarks

These timings come from `benchmarks/run.R`, run by the `benchmarks`
GitHub Actions workflow on one clean Ubuntu runner for the release and
for an earlier commit (`base`), one after the other. Times are medians
of three runs in seconds (one run for the slowest tasks); memory is what
one run allocates in R, measured with
[`bench::bench_memory()`](https://bench.r-lib.org/reference/bench_memory.html).
Solvers run with one thread. Run the script yourself with
`Rscript benchmarks/run.R`.

R 4.6.2, x86_64-pc-linux-gnu, optcompR 0.1.0.

## Building models

Construction paths are linear in the size of the model: concatenation,
indicator constraints, quadratic rows (stored as one triplet table) and
variables added one at a time no longer copy what is already there, and
a sparse coefficient matrix stays sparse in `A %*% x`.

| task | size | way | s (base) | s | speedup | MB (base) | MB |
|:---|:---|:---|---:|---:|---:|---:|---:|
| scalar constraints in a loop | 1000 | optcompR | 0.123 | 0.121 | 1.0 | 8.2 | 8.2 |
| one vector constraint | 1000 | optcompR | 0.001 | 0.002 | 0.5 | 0.5 | 0.5 |
| scalar constraints in a loop | 5000 | optcompR | 0.652 | 0.641 | 1.0 | 193.5 | 193.7 |
| one vector constraint | 5000 | optcompR | 0.003 | 0.004 | 0.8 | 2.5 | 2.7 |
| scalar variables in a loop | 2000 | optcompR | 0.176 | 0.155 | 1.1 | 162.1 | 3.0 |
| quadratic rows x_i^2 \<= 1 | 2000 | optcompR | 0.710 | 0.008 | 88.8 | 113.5 | 2.9 |
| indicators z_i =\> x_i \<= 1 | 2000 | optcompR | 0.175 | 0.016 | 10.9 | 233.3 | 3.6 |
| scalar variables in a loop | 8000 | optcompR | 2.001 | 0.657 | 3.0 | 2580.2 | 11.8 |
| quadratic rows x_i^2 \<= 1 | 8000 | optcompR | 4.349 | 0.028 | 155.3 | 1745.0 | 11.7 |
| indicators z_i =\> x_i \<= 1 | 8000 | optcompR | 4.588 | 0.061 | 75.2 | 3679.6 | 14.4 |

## Linear and mixed-integer programs

The same model solved through optcompR, through the solver’s own R
package and through ROI. The difference to the direct call is the cost
of building, classifying and checking the model.

| task | size | way | s (base) | s | speedup | MB (base) | MB |
|:---|:---|:---|---:|---:|---:|---:|---:|
| build | 1000 x 10000 | optcompR | 0.177 | 0.010 | 17.7 | 156.7 | 3.9 |
| build and solve | 1000 x 10000 | optcompR highs | 0.180 | 0.030 | 6.0 | 167.5 | 13.8 |
| build and solve | 1000 x 10000 | highs direct | 0.011 | 0.011 | 1.0 | 0.8 | 0.8 |
| build and solve | 1000 x 10000 | ROI highs | 0.017 | 0.017 | 1.0 | 2.7 | 2.7 |
| build | 5000 x 50000 | optcompR | 2.139 | 0.033 | 64.8 | 3835.3 | 19.5 |
| build and solve | 5000 x 50000 | optcompR highs | 2.240 | 0.166 | 13.5 | 3888.2 | 67.6 |
| build and solve | 5000 x 50000 | highs direct | 0.059 | 0.061 | 1.0 | 3.9 | 3.9 |
| build and solve | 5000 x 50000 | ROI highs | 0.075 | 0.075 | 1.0 | 13.4 | 13.4 |
| build and solve | 50 | optcompR highs | 0.037 | 0.036 | 1.0 | 0.1 | 0.1 |
| build and solve | 50 | optcompR scip | 0.015 | 0.014 | 1.1 | 0.1 | 0.1 |
| build and solve | 50 | highs direct | 0.031 | 0.029 | 1.1 | 0.0 | 0.0 |
| build and solve | 50 | ROI glpk | 0.005 | 0.005 | 1.0 | 0.1 | 0.1 |
| build and solve | 200 | optcompR highs | 0.016 | 0.015 | 1.1 | 0.3 | 0.3 |
| build and solve | 200 | optcompR scip | 0.012 | 0.011 | 1.1 | 0.3 | 0.3 |
| build and solve | 200 | highs direct | 0.010 | 0.009 | 1.1 | 0.0 | 0.0 |
| build and solve | 200 | ROI glpk | 0.010 | 0.009 | 1.1 | 0.2 | 0.2 |

## Quadratic and nonlinear programs

Automatic selection sends convex QPs to Clarabel: the HiGHS active set
QP solver needs thousands of iterations at a few thousand variables.

| task            | size | way               | s (base) |       s | speedup | MB (base) |   MB |
|:----------------|:-----|:------------------|---------:|--------:|--------:|----------:|-----:|
| build           | 1000 | optcompR          |    0.031 |   0.011 |     2.8 |      65.2 |  7.7 |
| solve           | 1000 | optcompR highs    |    0.427 |   0.396 |     1.1 |      21.2 | 21.4 |
| solve           | 1000 | optcompR osqp     |    0.085 |   0.058 |     1.5 |      21.2 | 21.5 |
| solve           | 1000 | optcompR clarabel |    0.115 |   0.114 |     1.0 |      21.3 | 21.6 |
| solve           | 1000 | osqp direct       |    0.010 |   0.010 |     1.0 |       0.3 |  0.3 |
| build           | 5000 | optcompR          |    1.229 |   0.039 |    31.5 |    1451.5 | 19.5 |
| solve           | 5000 | optcompR highs    |  101.174 | 100.649 |     1.0 |      27.7 | 29.0 |
| solve           | 5000 | optcompR osqp     |    1.752 |   0.909 |     1.9 |      27.8 | 29.1 |
| solve           | 5000 | optcompR clarabel |    5.606 |   5.591 |     1.0 |      28.3 | 29.6 |
| solve           | 5000 | osqp direct       |    0.473 |   0.473 |     1.0 |       1.3 |  1.3 |
| build and solve | 50   | optcompR ipopt    |    0.132 |   0.098 |     1.3 |       2.4 |  1.7 |
| build and solve | 50   | optcompR nloptr   |    0.312 |   0.277 |     1.1 |      15.0 | 12.5 |
| build and solve | 400  | optcompR ipopt    |    0.402 |   0.257 |     1.6 |      22.7 | 14.1 |

## Repeated solves and derivatives

After
[`opt_set_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_set_parameter.md),
lowered blocks whose parameters did not change are reused, and the HiGHS
solver and OSQP workspace of the previous solve are updated in place.
For the QP below, setting up a new OSQP workspace is cheap next to the
rest of each solve, so keeping it gains nothing measurable.

| task | size | way | s (base) | s | speedup | MB (base) | MB |
|:---|:---|:---|---:|---:|---:|---:|---:|
| 20 solves after parameter updates | 40 x 500 | optcompR parameters | 0.198 | 0.155 | 1.3 | 91.7 | 93.6 |
| 20 solves after parameter updates | 40 x 500 | rebuild each time | 0.371 | 0.348 | 1.1 | 133.5 | 128.2 |
| 20 scenarios | 40 x 500 | opt_scenarios | 0.197 | 0.154 | 1.3 | 91.7 | 93.6 |
| solution Jacobian | 20 | implicit KKT | 0.016 | 0.017 | 0.9 | 0.3 | 0.3 |
| solution Jacobian | 200 | implicit KKT | 0.023 | 0.023 | 1.0 | 10.9 | 8.0 |
| 20 QP solves after parameter updates | 2000 | optcompR osqp kept | 0.241 | 0.248 | 1.0 | 142.4 | 142.4 |
| 20 QP solves after parameter updates | 2000 | optcompR osqp new | 0.245 | 0.256 | 1.0 | 142.1 | 142.4 |
