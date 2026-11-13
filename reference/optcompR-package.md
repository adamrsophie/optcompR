# optcompR: Unified Mathematical Optimization and Computational Modeling in R

A solver independent modeling layer for mathematical optimization.
Decision variables, parameters, expressions, constraints and objectives
are written with ordinary R operators; models are classified (linear,
quadratic, mixed integer, nonlinear) with recorded evidence, lowered to
sparse standard forms, and solved by installed backends such as 'HiGHS'
through the 'highs' package, 'OSQP', 'NLopt', 'SCIP', and 'ROI' plugins.
Results carry standardized termination statuses, duals with a documented
sign convention, independent feasibility and Karush Kuhn Tucker
diagnostics, exact derivatives, sensitivity reports and reproducibility
metadata. The design follows 'JuMP' and 'MathOptInterface'
<https://jump.dev>.

## See also

Useful links:

- <https://github.com/adamrsophie/optcompR>

- <https://adamrsophie.github.io/optcompR/>

- Report bugs at <https://github.com/adamrsophie/optcompR/issues>

## Author

**Maintainer**: Adam R Sophie <adamrsophie@gmail.com> \[copyright
holder\]

Authors:

- Adam R Sophie <adamrsophie@gmail.com> \[copyright holder\]
