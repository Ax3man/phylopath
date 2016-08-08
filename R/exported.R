#' Compare causal models in a phylogenetic context.
#'
#' @param models A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{DAG}.
#' @param data A \code{data.frame} with data.
#' @param tree A phylogenetic tree of class \code{pylo}.
#' @param cor_fun A function that creates a \code{corStruct} object, typically
#'   one of the cor function from the \code{ape}, such as \code{corBrownian},
#'   \code{corPagel} etc.
#'
#' @return A table with relevant statistics for each model, including the C
#'   statistic, the associated p-values, information criterions, and model
#'   weigths.
#' @export
phylo_path <- function(models, data, tree, cor_fun = ape::corPagel) {
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

  d <- data.frame(model = names(models), k = k, q = q, C = C, p = p, CICc = IC)
  d <- d[order(d$CICc), ]
  d$delta_CICc <- d$CICc - d$CICc[1]
  d$l <- l(d$delta_CICc)
  d$w <- w(d$l)

  dplyr::mutate_if(d, is.numeric, round, digits = 3)
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