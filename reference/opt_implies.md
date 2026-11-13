# Logical constraints on binary variables

Exact linear formulations of logical relations between binary variables.
R's `&&`, `||` and `!` cannot build symbolic constraints, so these
functions are the explicit alternative (`1 - x` is the negation of a
binary `x`).

## Usage

``` r
opt_implies(model, a, b, name = NULL)

opt_and(model, z, x, name = NULL)

opt_or(model, z, x, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- a, b, z:

  Binary variables (scalars, or vectors of equal length for
  `opt_implies()`).

- name:

  Optional constraint name.

- x:

  A vector of binary variables.

## Value

An `opt_constraint` reference, invisibly.

## Details

- `opt_implies(model, a, b)`: \\a \Rightarrow b\\, i.e. \\a \le b\\; `a`
  and `b` may be negations such as `1 - x`.

- `opt_and(model, z, x)`: \\z = x_1 \wedge \dots \wedge x_k\\.

- `opt_or(model, z, x)`: \\z = x_1 \vee \dots \vee x_k\\.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 3, type = "binary")
z <- opt_variable(m, "z", type = "binary")
opt_and(m, z, x)
opt_constraint(m, sum(x) <= 2)
opt_objective(m, z)
opt_objective_value(opt_solve(m))  # 0: not all three can be 1
#> [1] 0
```
