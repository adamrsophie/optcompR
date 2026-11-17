# Read a CVXR solution back by variable name

Read a CVXR solution back by variable name

## Usage

``` r
opt_cvxr_values(cv)
```

## Arguments

- cv:

  An `opt_cvxr` object from
  [`opt_to_cvxr()`](https://adamrsophie.github.io/optcompR/reference/opt_to_cvxr.md)
  whose problem has been solved with
  [`CVXR::psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.html).

## Value

A named numeric vector, one value per model variable.
