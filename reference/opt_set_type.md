# Change variable types

Changing a type changes the problem class, for example from a linear to
a mixed-integer linear program. Setting `"binary"` also restricts the
bounds to `[0, 1]`.

## Usage

``` r
opt_set_type(model, x, type = c("continuous", "integer", "binary"))
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  A variable expression (or slice) or the variable block name.

- type:

  `"continuous"`, `"integer"` or `"binary"`.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 5)
opt_set_type(m, x, "integer")
opt_variables(m)
#>   name block    type lower upper start
#> 1    x     x integer     0     5    NA
```
