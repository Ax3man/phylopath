#' Compare causal models in a phylogenetic context.
#'
#' @param models A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{DAG}.
#' @param data A \code{data.frame} with data.
#' @param tree A phylogenetic tree of class \code{pylo}.
#' @param cor_fun A function that creates a \code{corStruct} object, typically
#'   one of the cor function from the \code{ape}, such as \code{corBrownian},
#'   \code{corPagel} etc.
#' @param cut_off The cut off value to decide which models to include in the
#'   average model, in terms of delta CICc.
#'
#' @return A table with relevant statistics for each model, including the C
#'   statistic, the associated p-values, information criterions, and model
#'   weigths.
#' @export
phylo_path <- function(models, data, tree, cor_fun = ape::corPagel,
                       cut_off = 2) {
  # Check if all models have the same number of nodes
  var_names <- lapply(models, colnames)
  if (stats::var(lengths(models)) != 0 |
      any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0)) {
    stop('All causal models need to include the same variables.', call. = FALSE)
  }

  if (is.null(names(models))) {
    names(models) <- LETTERS[1:length(models)]
  }
  formulas <- lapply(models, find_formulas)
  p_vals <- lapply(formulas, function(x) sapply(x, function(y) {
    m <- gls2(y, data = data, tree = tree, cor_fun = cor_fun)
    get_p(m)
  } ) )
  k <- lengths(formulas)
  q <- sapply(models, function(m) nrow(m) + sum(m))
  C <- sapply(p_vals, C_stat)
  p <- C_p(C, k)
  IC <- CICc(C, q, nrow(data))

  d <- data.frame(model = names(models), k = k, q = q, C = C, p = p, CICc = IC,
                  stringsAsFactors = FALSE)
  d <- d[order(d$CICc), ]
  d$delta_CICc <- d$CICc - d$CICc[1]
  d$l <- l(d$delta_CICc)
  d$w <- w(d$l)

  ps <- Map(function(x, y) {
    data.frame(d_sep = unlist(as.character(x)), p = unlist(y))
  }, formulas, p_vals)
  tab <- dplyr::mutate_if(d, is.numeric, round, digits = 3)

  best <- tab[tab$delta_CICc < 2, ]
  best_models <- lapply(models[best$model], est_DAG, data, cor_fun, tree)
  best_weigthed <- Map(`*`, best_models, best$w / sum(best$w))
  average <- apply(simplify2array(best_weigthed), c(1, 2), sum)
  class(average) <- c('matrix', 'DAG')

  out <- list(model_comp = tab,
              p_vals = ps,
              best_model = best_models[[1]],
              average_model = average)
  class(out) <- 'phylopath'
  return(out)
}

#' Directed acyclic graphs (DAGs)
#'
#' This function is a simple wrapper around the function from the \code{ggm}
#' package. The only differences are that the \code{order} argument defaults
#' to \code{TRUE} and that it adds a \code{DAG} class for easy plotting.
#'
#' @inheritParams ggm::DAG
#' @return An object of classes \code{matrix} and \code{DAG}
#' @export
DAG <- function(..., order = TRUE) {
  d <- ggm::DAG(..., order = order)
  class(d) <- c(class(d), 'DAG')
  d
}

#' Add normalized path coefficients to a DAG.
#'
#' @param DAG A directed acyclic graph, typically created with \code{DAG}.
#'
#' @inheritParams phylo_path
#' @return A matrix.
#' @export
est_DAG <- function(DAG, data, cor_fun, tree) {
  r <- rownames(data)
  data <- dplyr::mutate_if(data, is.numeric, scale)
  rownames(data) <- r
  d <- mapply(function(x, y, n) {
    if (all(y == 0)) {
      return(y)
    }
    f <- stats::formula(paste(x, paste(n[y == 1], collapse = '+'), sep = '~'))
    m <- gls2(f, data = data, cor_fun = cor_fun, tree = tree)
    y[y != 0] <- summary(m)$tTable[-1, 'Value']
    return(y)
  }, colnames(DAG), as.data.frame(DAG), MoreArgs = list(n = rownames(DAG)))
  class(d) <- c(class(d), 'DAG')
  d
}

#' @export
plot.DAG <- function(x, width_const = 5, ...) {
  df <- igraph::as_data_frame(
    igraph::graph_from_adjacency_matrix(x, weighted = TRUE), what = "both")
  df$vertices <- cbind(nodes = rownames(df$vertices),
                       df$vertices)
  if (!all(x == 0 | x == 1)) {
    df$edges$label <- round(df$edges$weight, 3)
    df$edges$penwidth <- abs(df$edges$weight / max(df$edges$weight) * width_const)
    df$edges$color <- ifelse(sign(df$edges$weight) == -1, 'red4', 'green4')
  }

  dg <- DiagrammeR::create_graph(
    nodes_df = df$vertices,
    edges_df = df$edges
  )
  DiagrammeR::render_graph(dg)
}

#' @export
print.phylopath <- function(x, ...) {
  plot.DAG(x$average_model)
  print(x$model_comp)
}

# function(x, node_df = NULL, edge_df = NULL) {
#   df <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(x, weighted = TRUE),
#                               what = "edges")
#   df$value <- df$weight
#   df$arrows <- 'to'
#   df$color <- 'black'
#
#   if (is.null(node_df)) {
#     node_df <- data.frame(id = rownames(x), label = rownames(x),
#                           shape = 'box', border = 'black', background = 'white')
#   }
#   if (!is.null(edge_df)) {
#     df <- merge(df, edge_df)
#   }
#   visNetwork::visNetwork(node_df, df)
# }
