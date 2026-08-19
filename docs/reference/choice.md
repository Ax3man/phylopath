# Extract and estimate an arbitrary model from a phylogenetic path analysis.

Extract and estimate an arbitrary model from a phylogenetic path
analysis.

## Usage

``` r
choice(phylopath, choice, ...)
```

## Arguments

- phylopath:

  An object of class `phylopath`.

- choice:

  A character string of the name of the model to be chosen, or the index
  in `model_set`.

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
    A = NL ~ BM,
    B = NL ~ LS,
    .common = c(LS ~ BM, DD ~ NL)
  )
  p <- phylo_path(candidates, rhino, rhino_tree)
  my_model <- choice(p, "B")
  # Print the best model to see coefficients, se and ci:
  my_model
#> A fitted causal model: 4 variables, 3 paths.
#>  Continuous:  BM LS NL DD 
#> 
#> Paths — standardized regression coefficients
#>     path coefficient    se
#>  BM → LS       0.497 0.089
#>  LS → NL       0.207 0.094
#>  NL → DD       0.629 0.080
  # Plot to show the weighted graph:
  plot(my_model)

```
