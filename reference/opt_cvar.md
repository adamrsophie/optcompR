# Risk measures and expectations over scenarios

- `opt_expectation(losses, prob)`: the sample average (or weighted
  expectation) of scenario losses, the objective of sample average
  approximation.

- `opt_cvar(model, losses, alpha, prob)`: conditional value at risk at
  level `alpha`, through the Rockafellar and Uryasev formulation \\t +
  \frac{1}{1-\alpha} \sum_k p_k (L_k - t)\_+\\ with auxiliary variables.
  Minimizing it, or bounding it from above, is exact.

- `opt_worst_expectation(model, losses, probs)`: distributionally robust
  expectation \\\max_j p_j^\top L\\ over a finite ambiguity set of
  distributions (rows of `probs`), and so over its convex hull, by an
  epigraph variable. Minimizing it, or bounding it from above, is exact.

## Usage

``` r
opt_cvar(model, losses, alpha, prob = NULL, name = NULL)

opt_expectation(losses, prob = NULL)

opt_worst_expectation(model, losses, probs, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- losses:

  Expression with one element per scenario.

- alpha:

  Confidence level in (0, 1).

- prob:

  Scenario probabilities (default uniform).

- name:

  Base name for generated variables.

- probs:

  Matrix of candidate distributions, one per row.

## Value

A scalar expression.

## Details

The second and third return expressions built on new variables and
constraints; they model the risk measure only in the direction stated
(in a minimization objective or on the left of `<=`).

## Examples

``` r
# portfolio over 4 return scenarios, minimizing CVaR of the loss
R <- matrix(c(0.10, 0.02, -0.05, 0.04,
              0.03, 0.01,  0.02, 0.01), 4)
m <- opt_model()
w <- opt_variable(m, "w", dim = 2, lower = 0)
opt_constraint(m, sum(w) == 1)
loss <- -(R %*% w)
opt_objective(m, opt_cvar(m, loss, alpha = 0.75))
opt_value(opt_solve(m), w)
#> [1] 0.125 0.875
```
