# Extract and average the best supported models from a phylogenetic path analysis.

Extract and average the best supported models from a phylogenetic path
analysis.

## Usage

``` r
average(phylopath, cut_off = 2, avg_method = "conditional", ...)
```

## Arguments

- phylopath:

  An object of class `phylopath`.

- cut_off:

  The CICc cut-off used to select the best models. Use `Inf` to average
  over all models. Use the
  [`best()`](https://ax3man.github.io/phylopath/reference/best.md)
  function to only use the top model, or
  [`choice()`](https://ax3man.github.io/phylopath/reference/choice.md)
  to select any single model.

- avg_method:

  Either `"full"` or `"conditional"`. The methods differ in how they
  deal with averaging a path coefficient where the path is absent in
  some of the models. The full method sets the coefficient (and the
  variance) for the missing paths to zero, meaning paths that are
  missing in some models will shrink towards zero. The conditional
  method only averages over models where the path appears, making it
  more sensitive to small effects. Following von Hardenberg &
  Gonzalez-Voyer 2013, conditional averaging is set as the default.

- ...:

  Arguments to pass to
  [phylolm::phylolm](https://rdrr.io/pkg/phylolm/man/phylolm.html) and
  [phylolm::phyloglm](https://rdrr.io/pkg/phylolm/man/phyloglm.html).
  Provide `boot = K` parameter to enable bootstrapping, where `K` is the
  number of bootstrap replicates. If you specified other options in the
  original
  [phylo_path](https://ax3man.github.io/phylopath/reference/phylo_path.md)
  call you don't need to specify them again.

## Value

An object of class `fitted_DAG`.

## See also

[`est_DAG()`](https://ax3man.github.io/phylopath/reference/est_DAG.md)
for what the coefficients mean.

## Examples

``` r
  candidates <- define_model_set(
    A = NL ~ RS,
    B = RS ~ NL + BM,
    .common = c(LS ~ BM, DD ~ NL, NL ~ BM)
  )
  p <- phylo_path(candidates, rhino, rhino_tree)
  summary(p)
#>   model k  q    C     p CICc delta_CICc     l     w
#> A     A 6  9 6.11 0.910 26.1       0.00 1.000 0.639
#> B     B 5 10 4.78 0.905 27.3       1.14 0.566 0.361

  # Models A and B have similar support, so we may decide to take
  # their average.

  avg_model <- average(p)
  # Print the average model to see coefficients, se and ci:
  avg_model
#> A fitted causal model: 5 variables, 6 paths.
#>  Continuous:  BM RS NL LS DD 
#> 
#> Paths — standardized regression coefficients
#>     path coefficient           95% CI
#>  BM → LS       0.497 [ 0.322,  0.673]
#>  BM → RS      -0.409 [-0.589, -0.229]
#>  BM → NL       0.450 [ 0.301,  0.599]
#>  RS → NL       0.528 [ 0.416,  0.640]
#>  NL → RS       0.874 [ 0.689,  1.059]
#>  NL → DD       0.629 [ 0.472,  0.785]

  if (FALSE) { # \dontrun{
  # Plot to show the weighted graph:
  plot(avg_model)

  # One can see that an averaged model is not necessarily a DAG itself.
  # This model actually has a path in two directions.

  # Note that coefficients that only occur in one of the models become much
  # smaller when we use full averaging:

  coef_plot(avg_model)
  coef_plot(average(p, method = 'full'))
  } # }
```
