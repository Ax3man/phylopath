# Plot a directed acyclic graph with path coefficients.

Plot a directed acyclic graph with path coefficients.

## Usage

``` r
# S3 method for class 'fitted_DAG'
plot(
  x,
  type = "width",
  labels = NULL,
  algorithm = "sugiyama",
  manual_layout = NULL,
  text_size = 6,
  box_x = 12,
  box_y = 8,
  edge_width = 1.25,
  curvature = 0.02,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  arrow = grid::arrow(type = "closed", 18, grid::unit(15, "points")),
  colors = c("firebrick", "navy"),
  show.legend = TRUE,
  width_const = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `fitted_DAG`.

- type:

  How to express the weight of the path. Either `"width"`, or `"color"`.

- labels:

  An optional set of labels to use for the nodes. This should be a named
  vector, of the form `c(var1 = "label1", var2 = "label2")`. If left at
  \`NULL“, the variable names of the DAGs are used.

- algorithm:

  A layout algorithm from `igraph`, see
  [`ggraph::create_layout()`](https://ggraph.data-imaginist.com/reference/ggraph.html).
  By default, uses the Sugiyama layout algorithm, which is designed to
  minimize edge crossing in DAGs.

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

- arrow:

  A [`grid::arrow`](https://rdrr.io/r/grid/arrow.html) object,
  specifying the shape and size of the arrowheads.

  The order of facets is taken from the ordering of the list, with the
  facet labels coming from the names of the list. If the list is
  unnamed, sequential lettering is used.

- colors:

  The end points of the continuous color scale. Keep in mind that red
  and green are obvious colors to use, but are better to be avoided
  because of color blind users.

- show.legend:

  Whether a legend for the color scale should be shown.

- width_const:

  Deprecated.

- ...:

  Not used.

## Examples

``` r
  d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
  d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
  plot(d_fitted)
```
