# Hessian of one quadratic row of a standard form

Hessian of one quadratic row of a standard form

## Usage

``` r
opt_row_hessian(sf, row)
```

## Arguments

- sf:

  An
  [`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md).

- row:

  Row index.

## Value

The symmetric sparse matrix \\Q_r\\ with quadratic part \\\frac12 x^\top
Q_r x\\, or `NULL` when the row has no quadratic terms.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
opt_constraint(m, x[1]^2 + x[1] * x[2] <= 1)
opt_row_hessian(opt_standard_form(m), 1)
#> 2 x 2 sparse Matrix of class "dsCMatrix"
#>         
#> [1,] 2 1
#> [2,] 1 .
```
