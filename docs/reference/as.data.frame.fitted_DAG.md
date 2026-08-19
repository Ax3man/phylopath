# Extract the paths of a fitted causal model.

Returns the estimated paths of a `fitted_DAG` as a long `data.frame`,
one row per path, which is a convenient starting point for custom tables
and figures. Paths that are absent from the causal model are not
included.

## Usage

``` r
# S3 method for class 'fitted_DAG'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., order_by = "default")
```

## Arguments

- x:

  An object of class `fitted_DAG`.

- row.names, optional, ...:

  Ignored, present for consistency with the generic.

- order_by:

  Either `"default"`, to order the paths by their position in the causal
  model, `"causal"`, to follow the paths downstream, or `"strength"`, to
  order them by the size of their coefficients. Ordering by strength is
  not available for models that contain both binary and continuous
  variables, since their coefficients are not comparable, see
  [`est_DAG()`](https://ax3man.github.io/phylopath/reference/est_DAG.md).
  Ordering by cause is not available for cyclical models, which
  [`average()`](https://ax3man.github.io/phylopath/reference/average.md)
  can produce when the direction of a path is not resolved.

## Value

A `data.frame` with columns `from`, `to`, `coef` and `se`, as well as
`lower` and `upper` if the model was fitted with `boot` larger than 0.

## Examples

``` r
  d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
  d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
  as.data.frame(d_fitted)
#>   from to        coef         se
#> 1   BM NL  0.43034428 0.08868688
#> 2   BM LS  0.49739369 0.08934185
#> 3   NL DD  0.63125366 0.08358289
#> 4   LS DD -0.01185413 0.08233878
  as.data.frame(d_fitted, order_by = 'strength')
#>   from to        coef         se
#> 1   LS DD -0.01185413 0.08233878
#> 2   BM NL  0.43034428 0.08868688
#> 3   BM LS  0.49739369 0.08934185
#> 4   NL DD  0.63125366 0.08358289
```
