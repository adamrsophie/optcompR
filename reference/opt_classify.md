# Classify an optimization model

Determines the problem class from the lowered standard form (after
parameter substitution) and records how convexity was established. Four
levels of evidence are distinguished: *syntactic* recognition of the
class, *structural* proofs (for example a diagonal quadratic form with
nonnegative entries), *numerical* evidence (eigenvalues compared with a
tolerance, which is not a formal proof near the boundary), and
properties that cannot be established automatically (general nonlinear
expressions are reported as `"unknown"`).

## Usage

``` r
opt_classify(model, tol = 1e-09)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md)
  or an `opt_standard_form`.

- tol:

  Relative eigenvalue tolerance for curvature checks.

## Value

An `opt_classification` list with `class`, `integer`, `objective`,
`constraints`, `convexity` (`"convex"`, `"nonconvex"` or `"unknown"`;
for mixed-integer models this refers to the continuous relaxation) and
`evidence`, a character vector.

## Details

Classes are `"LP"`, `"MILP"`, `"QP"`, `"MIQP"`, `"QCQP"`, `"MIQCQP"`,
`"NLP"` and `"MINLP"`. A model without an objective is a feasibility
problem and is classified by its constraints.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
opt_objective(m, x^2 + x * y + y^2)
opt_classify(m)
#> <opt_classification: QP, convex>
#>   - syntactic: objective is quadratic; constraint rows: none
#>   - numerical eigenvalues: objective is convex (smallest eigenvalue 1 (rounding level 1.33227e-13, tolerance 3e-09))
opt_objective(m, x^2 - 3 * x * y + y^2)
opt_classify(m)$convexity
#> [1] "nonconvex"
```
