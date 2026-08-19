# Plot the comparison of a set of causal models.

Shows the relative support for each causal model as a bar chart of model
weights, ordered from best to worst. Models within `cut_off` CICc of the
best model are highlighted, since those are usually considered to have
comparable support. Each bar is labelled with the p-value of the
d-separation test, where a significant value means the model is rejected
by the data.

## Usage

``` r
# S3 method for class 'phylopath_summary'
plot(x, cut_off = 2, ...)
```

## Arguments

- x:

  A `phylopath_summary` object, obtained by calling
  [`summary()`](https://rdrr.io/r/base/summary.html) on the result of
  [`phylo_path()`](https://ax3man.github.io/phylopath/reference/phylo_path.md).

- cut_off:

  The CICc difference within which models are highlighted as having
  support comparable to the best model.

- ...:

  Not used.

## Value

A `ggplot` object.

## Examples

``` r
  candidates <- define_model_set(
    A = NL ~ BM,
    B = NL ~ LS,
    .common = c(LS ~ BM, DD ~ NL)
  )
  p <- phylo_path(candidates, rhino, rhino_tree)
  plot(summary(p))
```
