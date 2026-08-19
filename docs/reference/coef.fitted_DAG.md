# Extract the coefficients and confidence intervals of a fitted causal model.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the path
coefficients as a named vector, and
[`confint()`](https://rdrr.io/r/stats/confint.html) the bounds of their
confidence intervals. Both name each path as `"from -> to"` and return
them in the same order, so they can be combined with
[`cbind()`](https://rdrr.io/r/base/cbind.html). Paths that are absent
from the causal model are not included. For a data frame with the
standard errors alongside the coefficients, see
[`as.data.frame.fitted_DAG()`](https://ax3man.github.io/phylopath/reference/as.data.frame.fitted_DAG.md).

## Usage

``` r
# S3 method for class 'fitted_DAG'
coef(object, ...)

# S3 method for class 'fitted_DAG'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `fitted_DAG`.

- ...:

  Ignored, present for consistency with the generic.

- parm:

  The paths for which to return intervals, given as names or indices.
  Defaults to all of them.

- level:

  The confidence level, which has to be 0.95. Only the bounds of the
  bootstrap interval are stored, not the replicates they were taken
  from, so another level cannot be produced after fitting.

## Value

For [`coef()`](https://rdrr.io/r/stats/coef.html), a named numeric
vector with one element per path. For
[`confint()`](https://rdrr.io/r/stats/confint.html), a matrix with a row
per path and the lower and upper bound in its two columns.

## Details

Note that these intervals are bootstrap percentile intervals, collected
while fitting, and not the parametric intervals that
[`confint()`](https://rdrr.io/r/stats/confint.html) returns for most
model objects.

## Examples

``` r
  d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
  d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
  coef(d_fitted)
#>    BM -> NL    BM -> LS    NL -> DD    LS -> DD 
#>  0.43034428  0.49739369  0.63125366 -0.01185413 

  # Confidence intervals require the model to be fitted with bootstrapping.
  # \donttest{
    d_boot <- est_DAG(d, rhino, rhino_tree, 'lambda', boot = 100)
    confint(d_boot)
#>               2.5 %    97.5 %
#> BM -> NL  0.2882040 0.5875534
#> BM -> LS  0.3105099 0.6676275
#> NL -> DD  0.4704635 0.7977846
#> LS -> DD -0.2205770 0.1150009
    cbind(coef = coef(d_boot), confint(d_boot))
#>                 coef      2.5 %    97.5 %
#> BM -> NL  0.43034428  0.2882040 0.5875534
#> BM -> LS  0.49739369  0.3105099 0.6676275
#> NL -> DD  0.63125366  0.4704635 0.7977846
#> LS -> DD -0.01185413 -0.2205770 0.1150009
  # }
```
