# Arithmetic and comparison operators for expressions

Expressions support `+`, `-`, `*`, `/` and `^` elementwise, with numeric
vectors, matrices and other expressions of the same model. Only
length-one operands are recycled; other length mismatches are errors.
The comparisons `<=`, `>=` and `==` build constraint relations for
[`opt_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_constraint.md);
they never return `TRUE` or `FALSE`. Strict inequalities `<` and `>` and
`!=` are rejected because they do not define closed feasible sets.
Logical operators are not overloaded; R's `&&` and `||` cannot be
overloaded at all, and logical constraints between expressions are not
supported.

## Usage

``` r
# S3 method for class 'opt_expr'
Ops(e1, e2)

# S3 method for class 'opt_expr'
Math(x, ...)

# S3 method for class 'opt_expr'
Summary(..., na.rm = FALSE)

# S3 method for class 'opt_expr'
matrixOps(x, y)

# S4 method for class 'Matrix,opt_expr'
x %*% y

# S4 method for class 'opt_expr,Matrix'
x %*% y
```

## Arguments

- e1, e2:

  Expressions or numeric values.

- x:

  An expression.

- ...:

  For [`log()`](https://rdrr.io/r/base/Log.html), an optional `base`;
  otherwise unused.

- na.rm:

  Ignored; expressions have no missing values.

- y:

  An expression or numeric matrix.

## Value

An `opt_expr` for arithmetic, an `opt_relation` for comparisons.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
2 * x + 1
#> <opt_expr 2, affine>
#>   [1] 2*x[1] + 1
#>   [2] 2*x[2] + 1
sum(x^2)
#> <opt_expr 1, quadratic>
#>   x[1]^2 + x[2]^2 
x[1] + x[2] <= 4
#> <opt_relation: 1 row, affine <= 0>
#>   x[1] + x[2] <= 4
```
