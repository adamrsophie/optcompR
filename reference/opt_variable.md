# Add decision variables to a model

Creates a scalar, vector, or array block of decision variables and
returns an expression that refers to them. Variables can be continuous,
integer or binary. Bounds are validated: contradictory bounds, integer
bounds without an integer between them, and binary bounds outside
`[0, 1]` are errors.

## Usage

``` r
opt_variable(
  model,
  name,
  dim = 1L,
  type = c("continuous", "integer", "binary"),
  lower = NULL,
  upper = NULL,
  dimnames = NULL,
  start = NULL
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- name:

  Block name. Must be unique among the model's variables and parameters.
  Element names are `name[i]` or `name[i,j]`, or use `dimnames`.

- dim:

  Integer dimensions: `1` for a scalar, `n` for a vector, `c(n, p)` for
  a matrix. Taken from `dimnames` when omitted.

- type:

  `"continuous"`, `"integer"` or `"binary"`.

- lower, upper:

  Bounds, recycled from length one. Binary variables have bounds within
  `[0, 1]`; the defaults for other types are free.

- dimnames:

  Optional character labels, a vector for one dimension or a list with
  one entry per dimension.

- start:

  Optional start values (used by solvers that accept warm starts).

## Value

An `opt_expr` with the shape given by `dim`.

## See also

[`opt_set_bounds()`](https://adamrsophie.github.io/optcompR/reference/opt_set_bounds.md),
[`opt_set_type()`](https://adamrsophie.github.io/optcompR/reference/opt_set_type.md),
[`opt_fix()`](https://adamrsophie.github.io/optcompR/reference/opt_fix.md)

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", dim = 3, type = "integer", lower = 0, upper = 10)
z <- opt_variable(m, "z", type = "binary",
                  dimnames = list(c("a", "b"), c("u", "v")))
opt_variables(m)
#>     name block       type lower upper start
#> 1      x     x continuous     0   Inf    NA
#> 2   y[1]     y    integer     0    10    NA
#> 3   y[2]     y    integer     0    10    NA
#> 4   y[3]     y    integer     0    10    NA
#> 5 z[a,u]     z     binary     0     1    NA
#> 6 z[b,u]     z     binary     0     1    NA
#> 7 z[a,v]     z     binary     0     1    NA
#> 8 z[b,v]     z     binary     0     1    NA
```
