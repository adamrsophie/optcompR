# Delete a constraint block

Delete a constraint block

## Usage

``` r
opt_delete_constraint(model, constraint)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- constraint:

  An `opt_constraint` reference or the block name.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
cn <- opt_constraint(m, x >= 1)
opt_delete_constraint(m, cn)
opt_constraints(m)
#> [1] name rows kind
#> <0 rows> (or 0-length row.names)
```
