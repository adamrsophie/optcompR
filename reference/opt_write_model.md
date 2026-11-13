# Write a model to an MPS or LP file

Writes the standard form of a linear, quadratic or mixed-integer linear
model in a format other solvers read, using the HiGHS writer. The format
follows the file extension: `.mps` for (free) MPS or `.lp` for the CPLEX
LP format. Parameters are written with their current values. Variable
and row names are not transferred; columns and rows keep the order of
[`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md).

## Usage

``` r
opt_write_model(model, file)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md)
  of class LP, MILP or QP.

- file:

  Path ending in `.mps` or `.lp`.

## Value

The file path, invisibly.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 3)
opt_constraint(m, x[1] + 2 * x[2] <= 4)
opt_objective(m, 3 * x[1] + x[2])
f <- tempfile(fileext = ".mps")
opt_write_model(m, f)
#> WARNING: Column names are not present, or contain duplicates: using names with prefix "c", beginning with suffix 0
#> WARNING: Row    names are not present, or contain duplicates: using names with prefix "r", beginning with suffix 0
#> Writing the model to /tmp/RtmpacUaT8/file1e7be96b2e7.mps
head(readLines(f))
#> [1] "NAME        " "OBJSENSE"     "  MAX"        "ROWS"         " N  Obj     "
#> [6] " L  r0      "
```
