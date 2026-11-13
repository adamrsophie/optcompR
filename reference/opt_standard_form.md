# Standard form of a model

Substitutes parameter values and lowers the model to sparse matrices in
the form \\\min / \max\\ \frac12 x^\top Q x + c^\top x +
\mathrm{offset}\\ subject to \\l_r \le A x + q_i(x) \le u_r\\, \\l \le x
\le u\\ and integrality. Rows keep a mapping back to the user's
constraint blocks. Results are cached and the cache is invalidated by
any change to the model structure, bounds or parameter values.

## Usage

``` r
opt_standard_form(model)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

## Value

An `opt_standard_form` list with elements `sense`, `c`, `Q`, `offset`,
`A`, `row_lower`, `row_upper`, `lower`, `upper`, `type`, `col_names`,
`row_names`, `row_block`, `row_kind` (`"linear"`, `"quadratic"` or
`"nonlinear"`), `row_Q` (list of row Hessian matrices in the one-half
convention for quadratic rows), and the lowered expressions used by
nonlinear solvers.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
opt_constraint(m, x + y <= 4)
opt_objective(m, 3 * x + 2 * y + 1)
sf <- opt_standard_form(m)
sf$c
#> [1] 3 2
as.matrix(sf$A)
#>      [,1] [,2]
#> [1,]    1    1
sf$offset
#> [1] 1
```
