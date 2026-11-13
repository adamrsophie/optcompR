# Symmetric matrix variables

Creates \\k(k+1)/2\\ decision variables and returns them as a symmetric
\\k \times k\\ matrix expression, for use with
[`opt_psd_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
and matrix products.

## Usage

``` r
opt_sym_variable(model, name, k, lower = NULL, upper = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- name:

  Block name. Must be unique among the model's variables and parameters.
  Element names are `name[i]` or `name[i,j]`, or use `dimnames`.

- k:

  Matrix size.

- lower, upper:

  Bounds, recycled from length one. Binary variables have bounds within
  `[0, 1]`; the defaults for other types are free.

## Value

A `k` by `k` `opt_expr`. The underlying variable block is named `name`
and holds the upper triangle by columns.

## Examples

``` r
# smallest t with t * I - A positive semidefinite: the largest eigenvalue
A <- matrix(c(2, 1, 1, 3), 2)
m <- opt_model()
X <- opt_sym_variable(m, "X", 2)
t <- opt_variable(m, "t")
opt_constraint(m, X == t * diag(2) - A)
opt_psd_constraint(m, X)
opt_objective(m, t)
opt_objective_value(opt_solve(m))
#> [1] 3.618034
max(eigen(A)$values)
#> [1] 3.618034
```
