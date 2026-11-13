# optcompR

Unified mathematical optimization and computational modeling in R.

optcompR is a solver independent modeling layer. Variables, parameters,
constraints and objectives are written with ordinary R operators; models are
classified (LP, QP, MILP, MIQP, QCQP, NLP, MINLP) with recorded evidence,
lowered to sparse standard forms, and solved by installed backends. Results
carry standardized statuses, duals with one documented sign convention, and
diagnostics recomputed independently from the model.

```r
library(optcompR)

m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
cap <- opt_constraint(m, x + y <= 4, name = "capacity")
opt_constraint(m, x <= 2)
opt_constraint(m, y <= 3)
opt_objective(m, 3 * x + 2 * y)

r <- opt_solve(m)
opt_status(r)          # "optimal"
opt_value(r, c(x, y))  # 2 2
opt_dual(r, cap)       # shadow price 2
opt_check(r)           # independent feasibility and KKT checks
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/optcompR")
```

## License

GPL-3.
