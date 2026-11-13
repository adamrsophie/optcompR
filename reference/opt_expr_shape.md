# Shape and indexing of expressions

Expressions behave like numeric vectors or matrices of scalar
expressions. Indexing accepts the same positive, negative, logical and
name indices as base R arrays. [`c()`](https://rdrr.io/r/base/c.html)
concatenates expressions of one model into a vector, and
[`t()`](https://rdrr.io/r/base/t.html) transposes a matrix expression.

## Usage

``` r
# S3 method for class 'opt_expr'
x[i, j, ..., drop = TRUE]

# S3 method for class 'opt_expr'
x[[i, ...]]

# S3 method for class 'opt_expr'
length(x)

# S3 method for class 'opt_expr'
dim(x)

# S3 method for class 'opt_expr'
dim(x) <- value

# S3 method for class 'opt_expr'
t(x)

# S3 method for class 'opt_expr'
c(...)
```

## Arguments

- x:

  An expression.

- i, j:

  Indices.

- ...:

  Further indices or expressions.

- drop:

  Drop dimensions of extent one, as in base R.

- value:

  New dimensions for `dim<-`.

## Value

An expression, or for [`length()`](https://rdrr.io/r/base/length.html)
and [`dim()`](https://rdrr.io/r/base/dim.html) integers.

## Examples

``` r
m <- opt_model()
X <- opt_variable(m, "X", dim = c(2, 3))
dim(X)
#> [1] 2 3
X[1, ]
#> <opt_expr 3, affine>
#>   [1] X[1,1]
#>   [2] X[1,2]
#>   [3] X[1,3]
t(X)
#> <opt_expr 3 x 2, affine>
#>   [1] X[1,1]
#>   [2] X[1,2]
#>   [3] X[1,3]
#>   [4] X[2,1]
#>   [5] X[2,2]
#>   [6] X[2,3]
c(X[1, 1], X[2, 3])
#> <opt_expr 2, affine>
#>   [1] X[1,1]
#>   [2] X[2,3]
```
