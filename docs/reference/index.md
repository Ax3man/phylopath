# Package index

## Top level functions

Use these functions to describe, fit, select and average your models.

- [`define_model_set()`](https://ax3man.github.io/phylopath/reference/define_model_set.md)
  : Define a model set.
- [`DAG()`](https://ax3man.github.io/phylopath/reference/DAG.md) :
  Directed acyclic graphs (DAGs)
- [`phylo_path()`](https://ax3man.github.io/phylopath/reference/phylo_path.md)
  : Compare causal models in a phylogenetic context.
- [`average()`](https://ax3man.github.io/phylopath/reference/average.md)
  : Extract and average the best supported models from a phylogenetic
  path analysis.
- [`best()`](https://ax3man.github.io/phylopath/reference/best.md) :
  Extract and estimate the best supported model from a phylogenetic path
  analysis.
- [`choice()`](https://ax3man.github.io/phylopath/reference/choice.md) :
  Extract and estimate an arbitrary model from a phylogenetic path
  analysis.
- [`show_warnings()`](https://ax3man.github.io/phylopath/reference/show_warnings.md)
  : Print out warnings from a phylopath analysis.

## Plotting

- [`plot(`*`<DAG>`*`)`](https://ax3man.github.io/phylopath/reference/plot.DAG.md)
  : Plot a directed acyclic graph.
- [`plot(`*`<fitted_DAG>`*`)`](https://ax3man.github.io/phylopath/reference/plot.fitted_DAG.md)
  : Plot a directed acyclic graph with path coefficients.
- [`plot_model_set()`](https://ax3man.github.io/phylopath/reference/plot_model_set.md)
  : Plot several causal hypothesis at once.
- [`coef_plot()`](https://ax3man.github.io/phylopath/reference/coef_plot.md)
  : Plot path coefficients and their confidence intervals or standard
  errors.

## Low level functions

Use these functions for finer control. These are mostly for internal
use.

- [`est_DAG()`](https://ax3man.github.io/phylopath/reference/est_DAG.md)
  : Add standardized path coefficients to a DAG.
- [`average_DAGs()`](https://ax3man.github.io/phylopath/reference/average_DAGs.md)
  : Perform model averaging on a list of DAGs.

## Included datasets

- [`rhino`](https://ax3man.github.io/phylopath/reference/rhino.md) :
  Rhinogrades traits.
- [`rhino_tree`](https://ax3man.github.io/phylopath/reference/rhino_tree.md)
  : Rhinogrades phylogeny.
- [`cichlids`](https://ax3man.github.io/phylopath/reference/cichlids.md)
  : Cichlid traits and the evolution of cooperative breeding.
- [`cichlids_tree`](https://ax3man.github.io/phylopath/reference/cichlids_tree.md)
  : Cichlid phylogeny.
- [`red_list`](https://ax3man.github.io/phylopath/reference/red_list.md)
  : Data on brain size, life history and vulnerability to extinction
- [`red_list_tree`](https://ax3man.github.io/phylopath/reference/red_list_tree.md)
  : Mammalian phylogeny
