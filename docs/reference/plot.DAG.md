# Plot a directed acyclic graph.

Plot a directed acyclic graph.

## Usage

``` r
# S3 method for class 'DAG'
plot(
  x,
  labels = NULL,
  algorithm = "sugiyama",
  manual_layout = NULL,
  text_size = 6,
  box_x = 12,
  box_y = 8,
  edge_width = 1.5,
  curvature = 0.02,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  arrow = grid::arrow(type = "closed", 18, grid::unit(15, "points")),
  ...
)
```

## Arguments

- x:

  A \`DAG“ object, usually created with the
  [`define_model_set()`](https://ax3man.github.io/phylopath/reference/define_model_set.md)
  or [`DAG()`](https://ax3man.github.io/phylopath/reference/DAG.md)
  function.

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

- ...:

  Not used.

## Examples

``` r
  d <- DAG(a ~ b + c + d)
  plot(d)


  # Plot with manually defined positions:
  ml <- data.frame(
    name = c('a', 'b', 'c', 'd'),
    x = c(1, 1, 2, 2),
    y = c(1, 2, 1, 2)
  )
  plot(d, manual_layout = ml)

```
