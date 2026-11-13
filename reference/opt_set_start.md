# Set start values for variables

Start values are passed to solvers that accept a warm start (for example
HiGHS for mixed-integer problems and NLopt for nonlinear problems).

## Usage

``` r
opt_set_start(model, x, value)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  A variable expression (or slice) or the variable block name.

- value:

  Start values, recycled from length one; `NA` clears them.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
opt_set_start(m, x, c(1, 2))
```
