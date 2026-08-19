#' @importFrom ggplot2 .data
NULL

#' @export
print.phylopath_summary <- function(x, ...) {
  print.data.frame(x, digits = 3)
  return(invisible(x))
}

# The paths of a fitted_DAG, formatted for printing.
path_table <- function(x, digits) {
  paths <- as.data.frame(x)
  fmt <- function(v) formatC(v, format = 'f', digits = digits)
  out <- data.frame(
    path = paste(paths$from, '\U2192', paths$to),
    coefficient = fmt(paths$coef)
  )
  if (!is.null(paths$lower)) {
    # Right aligning the bounds keeps the decimal points lined up, so a space
    # appears after the bracket only where another row needs room for a minus.
    out$`95% CI` <- paste0('[', format(fmt(paths$lower), justify = 'right'), ', ',
                           format(fmt(paths$upper), justify = 'right'), ']')
  } else {
    out$se <- fmt(paths$se)
  }
  out
}

#' @export
print.fitted_DAG <- function(x, digits = 3, ...) {
  n_paths <- sum(x$coef != 0)
  cat('A fitted causal model:', ncol(x$coef), 'variables,', n_paths,
      if (n_paths == 1) 'path.\n' else 'paths.\n')
  if (n_paths == 0) return(invisible(x))

  # Which variables are of which type, in the same terms print.phylopath() uses,
  # so that any single path can be interpreted. Empty categories are left out.
  if (!is.null(x$binary)) {
    used <- x$binary[used_vars(x$coef)]
    if (any(!used)) cat('\tContinuous:\t', names(used)[!used], '\n')
    if (any(used))  cat('\tBinary:\t\t', names(used)[used], '\n')
  }

  # Naming the scale is only possible when every coefficient is on the same one.
  # Otherwise the note below the table, and the warning, explain why not.
  scale <- coef_scale_name(x)
  cat('\n', if (is.null(scale)) 'Paths' else paste0('Paths \U2014 ', scale, 's'), '\n',
      sep = '')
  print.data.frame(path_table(x, digits), row.names = FALSE, right = TRUE)

  if (isTRUE(is_mixed(x))) {
    cat('\n', mixed_note(), '\n', sep = '')
    warn_if_mixed(x)
  }
  return(invisible(x))
}

#' Plot the comparison of a set of causal models.
#'
#' Shows the relative support for each causal model as a bar chart of model weights, ordered from
#' best to worst. Models within `cut_off` CICc of the best model are highlighted, since those are
#' usually considered to have comparable support. Each bar is labelled with the p-value of the
#' d-separation test, where a significant value means the model is rejected by the data.
#'
#' @param x A `phylopath_summary` object, obtained by calling `summary()` on the result of
#'   [phylo_path()].
#' @param cut_off The CICc difference within which models are highlighted as having support
#'   comparable to the best model.
#' @param ... Not used.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#'   candidates <- define_model_set(
#'     A = NL ~ BM,
#'     B = NL ~ LS,
#'     .common = c(LS ~ BM, DD ~ NL)
#'   )
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   plot(summary(p))
plot.phylopath_summary <- function(x, cut_off = 2, ...) {
  x$model <- factor(x$model, rev(x$model))
  ggplot2::ggplot(x, ggplot2::aes(x = .data$w, y = .data$model,
                                 fill = .data$delta_CICc < cut_off,
                                 label = round(.data$p, 3))) +
    ggplot2::geom_col(col = 'black', alpha = 0.6) +
    ggplot2::geom_text(hjust = "inward") +
    ggplot2::scale_fill_manual(
      values = c('TRUE' = 'firebrick', 'FALSE' = 'black'),
      labels = c('TRUE' = paste('within', cut_off, 'CICc')),
      breaks = c('TRUE')
    ) +
    ggplot2::scale_x_continuous(position = 'top') +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
    ggplot2::labs(y = "model weight", caption = "bar labels are p-values, significance indicates rejection") +
    ggplot2::theme(legend.position = 'bottom')
}

#' @export
print.phylopath <- function(x, ...) {
  num_vars <- names(x$data)[sapply(x$data, is.numeric)]
  bin_vars <- setdiff(names(x$data), num_vars)

  cat('A phylogenetic path analysis, on the variables:\n')
  # A category with no variables in it is left out, rather than printed empty.
  if (length(num_vars) > 0) cat('\tContinuous:\t', num_vars, '\n')
  if (length(bin_vars) > 0) cat('\tBinary:\t\t', bin_vars, '\n')
  cat('\n')
  cat(' Evaluated for these models:', names(x$model_set), '\n')
  cat('\n')
  cat(' Containing', sum(purrr::map_dbl(x$d_sep, nrow)), 'phylogenetic regressions, of which',
      length(unique(unlist(purrr::map(x$d_sep, 'd_sep')))), 'unique')
  cat('\n')
}

#' Print out warnings from a phylopath analysis.
#'
#' Use this function after running `phylo_path()` to conveniently print any generated warnings
#' to the screen. You can either provide no arguments, which will only work if you run it directly
#' after the analysis, or you have to provide the `phylopath` object manually.
#'
#' @param phylopath A `phylopath` object of which the warnings should be printed.
#'
#' @export
show_warnings <- function(phylopath = NULL) {
  if (is.null(phylopath)) phylopath <- .Last.value
  if (!inherits(phylopath, 'phylopath')) {
    rlang::abort(
      c('`phylopath` must be a phylopath object, as returned by `phylo_path()`.',
        x = paste0('It is of class ', paste(class(phylopath), collapse = '/'), '.'))
    )
  }
  sink <- lapply(phylopath$warnings, rlang::warn)
  return(invisible(NULL))
}

#' Plot a directed acyclic graph.
#'
#' @param x A `DAG`` object, usually created with the [define_model_set()] or [DAG()] function.
#' @param algorithm A layout algorithm from `igraph`, see
#'   [ggraph::create_layout()]. By default,
#'   uses the Sugiyama layout algorithm, which is designed to minimize edge crossing in DAGs.
#' @param ... Not used.
#' @inheritParams plot_model_set
#'
#' @export
#'
#' @examples
#'   d <- DAG(a ~ b + c + d)
#'   plot(d)
#'
#'   # Plot with manually defined positions:
#'   ml <- data.frame(
#'     name = c('a', 'b', 'c', 'd'),
#'     x = c(1, 1, 2, 2),
#'     y = c(1, 2, 1, 2)
#'   )
#'   plot(d, manual_layout = ml)
#'
plot.DAG <- function(x, labels = NULL, algorithm = 'sugiyama', manual_layout = NULL, text_size = 6,
                     box_x = 12, box_y = 8, edge_width = 1.5, curvature = 0.02, rotation = 0,
                     flip_x = FALSE, flip_y = FALSE,
                     arrow = grid::arrow(type = 'closed', 18, grid::unit(15, 'points')), ...) {
  if (sum(x) == 0) {
    rlang::abort('This causal model has no paths, so there is nothing to plot.')
  }

  g <- igraph::graph_from_adjacency_matrix(x, weighted = TRUE)

  l <- ggraph::create_layout(g, 'igraph', algorithm = algorithm)
  if (!is.null(manual_layout)) {
    l$x <- manual_layout$x[match(l$name, manual_layout$name)]
    l$y <- manual_layout$y[match(l$name, manual_layout$name)]
  }
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  ggraph::ggraph(l) +
    ggraph::geom_edge_arc(
      strength = curvature, arrow = arrow, edge_width = edge_width,
      end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
      start_cap = ggraph::rectangle(box_x, box_y, 'mm')
    ) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = text_size) +
    ggraph::theme_graph(base_family = 'sans')
}

#' Plot a directed acyclic graph with path coefficients.
#'
#' @param x An object of class `fitted_DAG`.
#' @param type How to express the weight of the path. Either `"width"`, or `"color"`.
#' @param algorithm A layout algorithm from \code{igraph}, see
#'   [ggraph::create_layout()]. By default,
#'   uses the Sugiyama layout algorithm, which is designed to minimize edge crossing in DAGs.
#' @param colors The end points of the continuous color scale. Keep in mind that red and green are
#'   obvious colors to use, but are better to be avoided because of color blind users.
#' @param show.legend Whether a legend for the color scale should be shown.
#' @param width_const Deprecated.
#' @param ... Not used.
#' @inheritParams plot_model_set
#'
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
#'   plot(d_fitted)
plot.fitted_DAG <- function(x, type = 'width', labels = NULL, algorithm = 'sugiyama',
                            manual_layout = NULL, text_size = 6, box_x = 12, box_y = 8,
                            edge_width = 1.25, curvature = 0.02, rotation = 0, flip_x = FALSE,
                            flip_y = FALSE,
                            arrow = grid::arrow(type = 'closed', 18, grid::unit(15, 'points')),
                            colors = c('firebrick', 'navy'), show.legend = TRUE,
                            width_const = NULL, ...) {
  if (!is.null(width_const)) {
    rlang::warn('The `width_const` argument is deprecated and is ignored.')
  }
  if (sum(x$coef) == 0) {
    rlang::abort('This causal model has no paths, so there is nothing to plot.')
  }
  type <- match.arg(type, c('width', 'color', 'colour'), FALSE)
  if (type == 'colour') type <- 'color'
  warn_if_mixed(x)

  # Name the scale rather than assuming it is standardized, and wrap it so that it
  # fits a legend key.
  scale <- coef_scale_name(x)
  scale_lab <- paste(strwrap(if (is.null(scale)) 'path coefficient' else scale,
                             width = 22), collapse = '\n')

  g <- igraph::graph_from_adjacency_matrix(x$coef, weighted = TRUE)
  l <- ggraph::create_layout(g, 'igraph', algorithm = algorithm)
  if (!is.null(manual_layout)) {
    l$x <- manual_layout$x[match(l$name, manual_layout$name)]
    l$y <- manual_layout$y[match(l$name, manual_layout$name)]
  }
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  if (type == 'width') {
    p <- ggplot2::ggplot(l) +
      ggraph::geom_edge_arc(
        ggplot2::aes(width = abs(.data$weight), color = .data$weight < 0,
                     label = round(.data$weight, 2)),
        strength = curvature, arrow = arrow, end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
        start_cap = ggraph::rectangle(box_x, box_y, 'mm'), show.legend = show.legend,
        linejoin = c('bevel'), angle_calc = 'along', label_dodge = grid::unit(10, 'points')) +
      ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = text_size) +
      ggraph::scale_edge_width_continuous(
        limits = c(0, max(abs(igraph::E(g)$weight))), range = c(0, 2), guide = 'none'
      ) +
      ggraph::scale_edge_color_manual(
        name = NULL,
        values = c('FALSE' = colors[2], 'TRUE' = colors[1]),
        labels = c('positive', 'negative')
      ) +
      ggraph::theme_graph(base_family = 'sans')
  }

  if (type == 'color') {
    p <- ggplot2::ggplot(l) +
      ggraph::geom_edge_arc(
        ggplot2::aes(colour = .data$weight, label = round(.data$weight, 2)),
        edge_width = edge_width,
        strength = curvature, arrow = arrow,
        end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
        start_cap = ggraph::rectangle(box_x, box_y, 'mm'),
        show.legend = show.legend,
        linejoin = c('bevel'),
        angle_calc = 'along',
        label_dodge = grid::unit(10, 'points')
      ) +
      ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = text_size) +
      ggraph::scale_edge_color_gradient2(
        scale_lab,
        low = colors[1], high = colors[2],
        limits = c(-max(abs(igraph::E(g)$weight)),
                   max(abs(igraph::E(g)$weight))),
        guide = ggraph::guide_edge_colorbar()
      ) +
      ggraph::theme_graph(base_family = 'sans')
  }
  if (isTRUE(is_mixed(x))) {
    p <- p + ggplot2::labs(caption = mixed_note())
  }
  return(p)
}

#' Plot path coefficients and their confidence intervals or standard errors.
#'
#' @param fitted_DAG A fitted DAG, usually obtained by [best()], [average()] or [est_DAG()].
#' @param error_bar Whether to use confidence intervals (`"ci"`) or standard errors (`"se"`) as
#'   error bars. Will force standard errors with a message if confidence intervals are not
#'   available.
#' @param order_by By `"default"`, the paths are ordered as in the the model that is supplied.
#'   Usually this is in the order that was established by `[phylo_path()]` for all combined graphs.
#'   This can be change to `"causal"` to do a reordering based on the model at hand, or to
#'   `"strength"` to order them by the standardized regression coefficient.
#' @param reverse_order If `TRUE`, the paths are plotted in reverse order.
#'   Particularly useful in combination with [ggplot2::coord_flip()] to create
#'   horizontal versions of the plot.
#' @param from Only show path coefficients from these nodes. Supply as a character vector.
#' @param to Only show path coefficients to these nodes. Supply as a character vector.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   plot(d)
#'   d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
#'   plot(d_fitted)
#'   coef_plot(d_fitted, error_bar = "se")
#'   # to create a horizontal version, use this:
#'   coef_plot(d_fitted, error_bar = "se", reverse_order = TRUE) + ggplot2::coord_flip()
coef_plot <- function(fitted_DAG, error_bar = 'ci', order_by = "default", from = NULL, to = NULL,
                      reverse_order = FALSE) {
  stopifnot(inherits(fitted_DAG, 'fitted_DAG'))
  error_bar <- match.arg(error_bar, c('ci', 'se'), several.ok = FALSE)
  order_by <- match.arg(order_by, c('default', 'causal', 'strength'), FALSE)
  warn_if_mixed(fitted_DAG)
  if (error_bar == 'ci' & is.null(fitted_DAG$lower)) {
    rlang::inform(
      c('This model has no confidence intervals, so standard errors are shown instead.',
        i = 'Fit the model with `boot` larger than 0 for intervals, or set `error_bar = "se"` to silence this message.')
    )
    error_bar <- 'se'
  }

  df <- as.data.frame(fitted_DAG, order_by = order_by)
  if (error_bar == 'se') {
    df$lower <- df$coef - df$se
    df$upper <- df$coef + df$se
  }
  df$path <- paste(df$from, df$to, sep = ' \U2192 ')

  # Do the filtering of paths:
  if (!is.null(from)) {
    df <- df[df$from %in% from, ]
  }
  if (!is.null(to)) {
    # make sure the variable is unambiguous
    df <- df[df$to %in% to, ]
  }

  if (reverse_order) {
    df$path <- factor(df$path, levels = rev(df$path))
  } else {
    df$path <- factor(df$path, levels = df$path)
  }

  # Name the scale on the axis when every coefficient is on the same one, which is
  # the case unless the model mixes binary and continuous variables.
  scale <- coef_scale_name(fitted_DAG)
  y_lab <- paste(if (is.null(scale)) 'path coefficient' else scale,
                 if (error_bar == 'ci') '\U00B1 CI' else '\U00B1 SE')

  p <- ggplot2::ggplot(df,
                  ggplot2::aes(x = .data$path, y = .data$coef,
                               ymin = .data$lower, ymax = .data$upper)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 1, lty = 2) +
    ggplot2::geom_pointrange(size = 0.75) +
    ggplot2::xlab('') +
    ggplot2::ylab(y_lab)
  if (isTRUE(is_mixed(fitted_DAG))) {
    p <- p + ggplot2::labs(caption = mixed_note())
  }
  p
}


#' Plot several causal hypothesis at once.
#'
#' @param model_set A list of `DAG` objects, usually created with [define_model_set()].
#' @param labels An optional set of labels to use for the nodes. This should be a named vector, of
#'   the form `c(var1 = "label1", var2 = "label2")`.
#'   If left at `NULL``, the variable names of the DAGs are used.
#' @param algorithm A layout algorithm from `igraph`, see
#'   [ggraph::create_layout()]. By default, uses the Kamada-Kawai
#'   layout algorithm. Another good option is `"sugiyama"`, which is
#'   designed to minimize edge crossing in DAGs. However, it can often plot
#'   nodes too close together.
#' @param manual_layout Alternatively, precisely define the layout yourself, by providing a
#'   `data.frame` that at least has a column `name` with all variable names, and columns `x` and `y`
#'   with positions to be plotted. Setting this parameter overrides `algorithm` but other changes,
#'   such as `rotation` and `flip`s will still be applied.
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
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' m <- list(one = DAG(a ~ b + c + d), two = DAG(a ~ b, b ~ c, d ~ d))
#' plot_model_set(m)
#' plot_model_set(m, algorithm = "sugiyama")
plot_model_set <- function(model_set, labels = NULL, algorithm = 'kk', manual_layout = NULL,
                           text_size = 5, box_x = 12, box_y = 10, edge_width = 1, curvature = 0.05,
                           rotation = 0, flip_x = FALSE, flip_y = FALSE, nrow = NULL,
                           arrow = grid::arrow(type = 'closed', 15, grid::unit(10, 'points'))) {
  # Input checks
  if (!is.list(model_set) | !all(purrr::map_lgl(model_set, ~inherits(., 'DAG')))) {
    rlang::abort(
      c('`model_set` must be a list of DAG objects.',
        i = 'Create one with `define_model_set()`.')
    )
  }
  names(model_set) <- name_model_set(names(model_set), length(model_set))
  stop_if_variables_differ(model_set)

  # Build  single complete graph
  result <- igraph::make_empty_graph() + igraph::vertices(row.names(model_set[[1]]))
  for (i in seq_along(model_set)) {
    m <- model_set[[i]]
    ind <- which(m == 1, arr.ind = TRUE)
    result <- igraph::add_edges(
      result,
      c(rbind(rownames(m)[ind[, 'row']], colnames(m)[ind[, 'col']])),
      attr = list(model = names(model_set)[[i]])
    )
  }
  igraph::edge.attributes(result)$model <- factor(igraph::E(result)$model,
                                                  names(model_set), names(model_set))

  l <- ggraph::create_layout(result, 'igraph', algorithm = algorithm)
  if (!is.null(manual_layout)) {
    l$x <- manual_layout$x[match(l$name, manual_layout$name)]
    l$y <- manual_layout$y[match(l$name, manual_layout$name)]
  }
  l <- adjust_layout(l, rotation, flip_x, flip_y)
  l <- combine_with_labels(l, labels)

  # Build plot.
  ggraph::ggraph(l) +
    ggraph::geom_edge_arc(strength = curvature, arrow = arrow, edge_width = edge_width,
                           end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
                           start_cap = ggraph::rectangle(box_x, box_y, 'mm')) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = text_size) +
    ggraph::facet_edges(~model, nrow = nrow) +
    ggplot2::scale_x_continuous(expand = c(0.2, 0)) +
    ggplot2::scale_y_continuous(expand = c(0.2, 0)) +
    ggraph::theme_graph(foreground = 'grey80', base_family = 'sans')
}
