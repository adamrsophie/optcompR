# Create an optimization model

An `opt_model` is a mutable container for decision variables,
parameters, constraints and one objective. Models have reference
semantics: functions such as
[`opt_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_variable.md)
and
[`opt_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_constraint.md)
change the model in place, and every copy of the object refers to the
same model. Expressions remember the model they were built from, and
combining expressions from two models is an error.

## Usage

``` r
opt_model(sense = c("min", "max"), name = NULL)
```

## Arguments

- sense:

  Optimization sense, `"min"` or `"max"`.

- name:

  Optional model name used in printing and reports.

## Value

An object of class `opt_model`.

## Lifetime and serialization

A model lives as long as some object refers to it. Saving a model with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) stores its variables,
constraints and objective. Expressions saved in a separate file refer to
their own copy of the model after loading, so save a model together with
the expressions you want to reuse, for example in one list.

## See also

[`opt_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_variable.md),
[`opt_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_constraint.md),
[`opt_objective()`](https://adamrsophie.github.io/optcompR/reference/opt_objective.md),
[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0, upper = 2)
y <- opt_variable(m, "y", lower = 0, upper = 3)
opt_constraint(m, x + y <= 4, name = "capacity")
opt_objective(m, 3 * x + 2 * y)
m
#> <opt_model>
#>   sense:       max
#>   variables:   2 (0 integer or binary)
#>   constraints: 1 rows in 1 blocks
#>   objective:   affine
```
