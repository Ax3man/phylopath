# Plot several causal hypothesis at once.

Plot several causal hypothesis at once.

## Usage

``` r
plot_model_set(
  model_set,
  labels = NULL,
  algorithm = "kk",
  manual_layout = NULL,
  text_size = 5,
  box_x = 12,
  box_y = 10,
  edge_width = 1,
  curvature = 0.05,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  nrow = NULL,
  arrow = grid::arrow(type = "closed", 15, grid::unit(10, "points"))
)
```

## Arguments

- model_set:

  A list of `DAG` objects, usually created with
  [`define_model_set()`](https://ax3man.github.io/phylopath/reference/define_model_set.md).

- labels:

  An optional set of labels to use for the nodes. This should be a named
  vector, of the form `c(var1 = "label1", var2 = "label2")`. If left at
  \`NULL“, the variable names of the DAGs are used.

- algorithm:

  A layout algorithm from `igraph`, see
  [`ggraph::create_layout()`](https://ggraph.data-imaginist.com/reference/ggraph.html).
  By default, uses the Kamada-Kawai layout algorithm. Another good
  option is `"sugiyama"`, which is designed to minimize edge crossing in
  DAGs. However, it can often plot nodes too close together.

- manual_layout:

  Alternatively, precisely define the layout yourself, by providing a
  `data.frame` that at least has a column `name` with all variable
  names, and columns `x` and `y` with positions to be plotted. Setting
  this parameter overrides `algorithm` but other changes, such as
  `rotation` and `flip`s will still be applied.

- text_size:

  Size of the node label text.

- box_x:

  To avoid the arrows colliding with the nodes, specify the rectangular
  dimensions of an invisible box around each node. If you have long
  labels, you need to increase this.

- box_y:

  To avoid the arrows colliding with the nodes, specify the rectangular
  dimensions of an invisible box around each node. If you have
  multi-line labels, you need to increase this.

- edge_width:

  Width of the edges.

- curvature:

  Curvature of the edges. A slight curvature can look pretty.

- rotation:

  Supply the degrees you want to rotate the layout by. This is useful in
  order to put rotate your upstream nodes towards the top if needed.

- flip_x:

  Whether to flip the node positions horizontally.

- flip_y:

  Whether to flip the node positions vertically.

- nrow:

  Number of rows to display the models on.

- arrow:

  A [`grid::arrow`](https://rdrr.io/r/grid/arrow.html) object,
  specifying the shape and size of the arrowheads.

  The order of facets is taken from the ordering of the list, with the
  facet labels coming from the names of the list. If the list is
  unnamed, sequential lettering is used.

## Value

A `ggplot` object.

## Examples

``` r
m <- list(one = DAG(a ~ b + c + d), two = DAG(a ~ b, b ~ c, d ~ d))
plot_model_set(m)

plot_model_set(m, algorithm = "sugiyama")
```
