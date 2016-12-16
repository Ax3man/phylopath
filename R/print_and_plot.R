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

#' @export
plot.DAG <- function(x, ...) {
  df <- igraph::as_data_frame(
    igraph::graph_from_adjacency_matrix(x, weighted = TRUE), what = "both")

  nodes_df <- DiagrammeR::create_node_df(n = nrow(df$vertices),
                                         label = df$vertices$name,
                                         shape = 'oval')
  edges_df <- DiagrammeR::create_edge_df(match(df$edges$from, df$vertices$name),
                                         match(df$edges$to, df$vertices$name))

  dg <- DiagrammeR::create_graph(
    nodes_df = nodes_df,
    edges_df = edges_df
  )
  DiagrammeR::render_graph(dg, title = '', ...)
}

#' @export
plot.fitted_DAG <- function(x, width_const = 5, ...) {
  df <- igraph::as_data_frame(
    igraph::graph_from_adjacency_matrix(x$coef, weighted = TRUE), what = "both")

  nodes_df <- DiagrammeR::create_node_df(n = nrow(df$vertices),
                                         label = df$vertices$name,
                                         shape = 'oval')
  edges_df <- DiagrammeR::create_edge_df(
    from = match(df$edges$from, df$vertices$name),
    to = match(df$edges$to, df$vertices$name),
    rel = round(df$edges$weight, 3),
    label = round(df$edges$weight, 3),
    penwidth = abs(df$edges$weight / max(df$edges$weight) * width_const),
    color = ifelse(sign(df$edges$weight) == -1, 'red4', 'green4')
  )

  dg <- DiagrammeR::create_graph(
    nodes_df = nodes_df,
    edges_df = edges_df
  )
  DiagrammeR::render_graph(dg, title = '', ...)
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
  df$path <- paste(df$from, df$to, sep = ' -> ')
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