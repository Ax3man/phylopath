#' @export
print.phylopath_summary <- function(x, ...) {
  print(dplyr::mutate_if(x, is.numeric, dplyr::funs(round), digits = 3))
  return(invisible(x))
}

#' @export
print.phylopath <- function(x, ...) {
  cat('\n')
  cat('A phylogenetic path analysis\n')
  cat('\n')
  cat('  Evaluated for these models:', names(x$models), '\n')
  cat('\n')
  cat('  Containing', sum(purrr::map_dbl(x$d_sep, nrow)), 'phylogenetic regressions.')
  cat('\n')
}

#' Plot a directed acyclic graph.
#'
#' @param x A \code{DAG} object, usually created with the \code{DAG} function.
#' @param algorithm A layout algorithm from \code{igraph}, see
#'   \code{\link[ggraph]{create_layout}} and \code{\link[ggraph]{create_layout.igraph}}. By default,
#'   uses the Sugiyama layout algorithm, which is designed to minimize edge crossing in DAGs.
#' @param ... Not used.
#' @inheritParams plot_model_set
#'
#' @export
plot.DAG <- function(x, labels = NULL, algorithm = 'sugiyama', text_size = 6, box_x = 12, box_y = 8,
                     edge_width = 1.5, curvature = 0.02, rotation = 0, flip_x = FALSE,
                     flip_y = FALSE,
                     arrow = grid::arrow(type = 'closed', 18, grid::unit(15, 'points')), ...) {
  g <- igraph::graph_from_adjacency_matrix(x, weighted = TRUE)

  l <- ggraph::create_layout(g, 'igraph', algorithm = algorithm)
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  ggraph::ggraph(l) +
    ggraph::geom_edge_arc(curvature = curvature, arrow = arrow, edge_width = edge_width,
                          end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
                          start_cap = ggraph::rectangle(box_x, box_y, 'mm')) +
    ggraph::geom_node_text(ggplot2::aes_(label = ~name), size = text_size) +
    ggraph::theme_graph(base_family = 'sans')
}

#' Plot a directed acyclic graph with path coefficients.
#'
#' @param x An object of class \code{fitted_DAG}.
#' @param algorithm A layout algorithm from \code{igraph}, see
#'   \code{\link[ggraph]{create_layout}} and \code{\link[ggraph]{create_layout.igraph}}. By default,
#'   uses the Sugiyama layout algorithm, which is designed to minimize edge crossing in DAGs.
#' @param colors The end points of the continious color scale. Keep in mind that red and green are
#'   obvious colors to use, but are better to be avoided because of color blind users.
#' @param show.legend Whether a legend for the color scale should be shown.
#' @param width_const Deprecated.
#' @param ... Not used.
#' @inheritParams plot_model_set
#'
#' @export
plot.fitted_DAG <- function(x, labels = NULL, algorithm = 'sugiyama', text_size = 6, box_x = 12,
                            box_y = 8, edge_width = 1.25, curvature = 0.02, rotation = 0,
                            flip_x = FALSE, flip_y = FALSE,
                            arrow = grid::arrow(type = 'closed', 18, grid::unit(15, 'points')),
                            colors = c('firebrick', 'navy'), show.legend = TRUE,
                            width_const = NULL, ...) {
  if (!is.null(width_const)) {
    warning('width_const has been deprecated and is ignored.', call. = FALSE)
  }

  g <- igraph::graph_from_adjacency_matrix(x$coef, weighted = TRUE)
  l <- ggraph::create_layout(g, 'igraph', algorithm = algorithm)
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  suppressWarnings( {
    ggplot2::ggplot(l) +
      ggraph::geom_edge_arc(ggplot2::aes_(colour = ~weight, label = ~round(weight, 2)),
                            edge_width = edge_width,
                            curvature = curvature, arrow = arrow,
                            end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
                            start_cap = ggraph::rectangle(box_x, box_y, 'mm'),
                            show.legend = show.legend,
                            linejoin = c('bevel'),
                            angle_calc = 'along',
                            label_dodge = grid::unit(10, 'points')) +
      ggraph::geom_node_text(ggplot2::aes_(label = ~name), size = text_size) +
      ggraph::scale_edge_color_gradient2('standardized\npath coefficient',
                                         low = colors[1], high = colors[2],
                                         limits = c(-max(abs(igraph::E(g)$weight)),
                                                    max(abs(igraph::E(g)$weight))),
                                         guide = ggraph::guide_edge_colorbar()) +
      ggraph::theme_graph(base_family = 'sans')
  } )
}

#' Plot path coefficients and their confidence intervals.
#'
#' @param fitted_DAG A fitted DAG, usually obtained by \code{best},
#'   \code{average} or \code{est_DAG}.
#' @param reverse_order If \code{TRUE}, the paths are plotted in reverse order.
#'   Particularly useful in combination with \code{ggplot2::coor_flip} to create
#'   horizontal versions of the plot.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   plot(d)
#'   d_fitted <- est_DAG(d, rhino, ape::corBrownian, rhino_tree)
#'   plot(d_fitted)
#'   coef_plot(d_fitted)
#'   # to create a horizontal version, use this:
#'   coef_plot(d_fitted, reverse_order = TRUE) + ggplot2::coord_flip()
coef_plot <- function(fitted_DAG, reverse_order = FALSE) {
  df <- as.data.frame(fitted_DAG$coef)
  df <- tibble::rownames_to_column(df, 'from')
  df <- tidyr::gather_(df, 'to', 'coef', colnames(fitted_DAG$coef))
  df$lower <- c(fitted_DAG$lower)
  df$upper <- c(fitted_DAG$upper)
  df$path <- paste(df$from, df$to, sep = ' \U2192 ')
  df <- dplyr::arrange_(df,
                        ~match(df$from, colnames(fitted_DAG$coef)),
                        ~match(df$to, colnames(fitted_DAG$coef)))
  if (reverse_order) {
    df$path <- factor(df$path, levels = rev(df$path))
  } else {
    df$path <- factor(df$path, levels = df$path)
  }
  df <- df[abs(df$coef) > .Machine$double.eps, ]
  ggplot2::ggplot(df,
                  ggplot2::aes_(~path, ~coef, ymin = ~lower, ymax = ~upper)) +
    ggplot2::geom_hline(yintercept = 0, size = 1, lty = 2) +
    ggplot2::geom_pointrange(size = 0.75) +
    ggplot2::xlab('') +
    ggplot2::ylab('standardized regression coefficient \U00B1 CI')
}


#' Plot several causal hypothesis at once.
#'
#' @param models A list of \code{DAG} objects.
#' @param labels An optional set of labels to use for the nodes. This should be a named vector, of
#'   the form \code{c(var1 = "label1", var2 = "label2")}.
#'   If left at \code{NULL}, the variable names of the DAGs are used.
#' @param algorithm A layout algorithm from \code{igraph}, see
#'   \code{\link[ggraph]{create_layout}}. By default, uses the Kamada-Kawai
#'   layout algorithm. Another good option is \code{"sugiyama"}, which is
#'   designed to minimize edge crossing in DAGs. However, it can often plot
#'   nodes too close together.
#' @param text_size Size of the node label text.
#' @param box_x To avoid the arrows colliding with the nodes, specify the
#'   rectangular dimensions of an invisible box around each node. If you have
#'   long labels, you need to increase this.
#' @param box_y To avoid the arrows colliding with the nodes, specify the
#'   rectangular dimensions of an invisible box around each node. If you have
#'   multi-line labels, you need to increase this.
#' @param edge_width Width of the edges.
#' @param curvature Curvature of the edges. A slight curvature can look pretty.
#' @param rotation Supply the degrees you want to rotate the layout by. This is useful in order to
#'   put rotate your upstream nodes towards the top if needed.
#' @param flip_x Whether to flip the node positions horizontally.
#' @param flip_y Whether to flip the node positions vertically.
#' @param nrow Number of rows to display the models on.
#' @param arrow A \code{grid::arrow} object, specifying the shape and size of the arrowheads.
#'
#' The order of facets is taken from the ordering of the list, with the facet
#' labels coming from the names of the list. If the list is unnamed, sequential
#' lettering is used.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' m <- list(one = DAG(a ~ b + c + d), two = DAG(a ~ b, b ~ c, d ~ d))
#' plot_model_set(m)
#' plot_model_set(m, algorithm = "sugiyama")
plot_model_set <- function(models, labels = NULL, algorithm = 'kk', text_size = 5, box_x = 12,
                           box_y = 10, edge_width = 1, curvature = 0.05, rotation = 0,
                           flip_x = FALSE, flip_y = FALSE, nrow = NULL,
                           arrow = grid::arrow(type = 'closed', 15, grid::unit(10, 'points'))) {
  # Input checks
  if (!is.list(models) | !all(purrr::map_lgl(models, ~inherits(., 'DAG')))) {
    stop('models should be a list of DAG objects.')
  }
  if (is.null(names(models))) {
    names(models) <- LETTERS[seq_along(models)]
  }
  var_names <- lapply(models, colnames)
  if (length(models) > 1 &
      (stats::var(lengths(models)) != 0 |
       any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0))) {
    stop('All causal models need to include the same variables. Combined, your
         models include the following variables:\n',
         paste(sort(unique(unlist(var_names))), collapse = '\n'),
         call. = FALSE)
  }

  # Build  single complete graph
  result <- igraph::make_empty_graph() + igraph::vertices(row.names(models[[1]]))
  for (i in seq_along(models)) {
    m <- models[[i]]
    ind  <- which(m == 1)
    from <- ind %% nrow(models[[i]])
    to   <- (ind - from) / nrow(models[[i]]) + 1
    result <- igraph::add_edges(result, c(rbind(rownames(m)[from], colnames(m)[to])),
                                attr = list(model = names(models)[[i]]))
  }
  igraph::edge.attributes(result)$model <- factor(igraph::E(result)$model,
                                                  names(models), names(models))

  l <- ggraph::create_layout(result, 'igraph', algorithm = algorithm)
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  # Build plot.
  ggraph::ggraph(l) +
    ggraph::geom_edge_arc(curvature = curvature, arrow = arrow, edge_width = edge_width,
                           end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
                           start_cap = ggraph::rectangle(box_x, box_y, 'mm')) +
    ggraph::geom_node_text(ggplot2::aes_(label = ~name), size = text_size) +
    ggraph::facet_edges(~model, nrow = nrow) +
    ggplot2::scale_x_continuous(expand = c(0.2, 0)) +
    ggplot2::scale_y_continuous(expand = c(0.2, 0)) +
    ggraph::theme_graph(foreground = 'grey80', base_family = 'sans')
}
