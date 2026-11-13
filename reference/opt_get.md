# Retrieve a variable or parameter block by name

Retrieve a variable or parameter block by name

## Usage

``` r
opt_get(model, name)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- name:

  Block name.

## Value

The expression for the whole block.

## Examples

``` r
m <- opt_model()
opt_variable(m, "x", dim = 3)
#> <opt_expr 3, affine>
#>   [1] x[1]
#>   [2] x[2]
#>   [3] x[3]
opt_get(m, "x")
#> <opt_expr 3, affine>
#>   [1] x[1]
#>   [2] x[2]
#>   [3] x[3]
```
