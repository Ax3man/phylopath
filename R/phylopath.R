#' Compare causal models in a phylogenetic context.
#'
#' @param models A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{DAG}.
#' @param data A \code{data.frame} with data.
#' @param tree A phylogenetic tree of class \code{pylo}.
#' @param cor_fun A function that creates a \code{corStruct} object, typically
#'   one of the cor function from the \code{ape}, such as \code{corBrownian},
#'   \code{corPagel} etc.
#' @param order Causal order of the included variable, given as a character
#'   vector. This is used to determine which variable should be the dependent
#'   in the dsep regression equations. If left unspecified, the order will be
#'   automatically determined. If the combination of all included models is
#'   itself a DAG, then the ordering of that full model is used. Otherwise,
#'   the most common ordering between each pair of variables is used to create
#'   a general ordering.
#' @param parallel An optional vector containing the virtual connection
#'   process type for running the chains in parallel (such as \code{"SOCK"}).
#'   A cluster is create using the \code{parallel} package.
#'
#' @return A phylopath object, with the following components:
#'  \describe{
#'   \item{d_sep}{for each model a table with seperation statements and statistics.}
#'   \item{models}{the DAGs}
#'   \item{data}{the supplied data}
#'   \item{tree}{the supplied tree}
#'   }
#' @export
#' @examples
#'   #see vignette('intro_to_phylopath') for more details
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'
#'   # Printing p gives some general information:
#'   p
#'   # And the summary gives statistics to compare the models:
#'   summary(p)
#'
phylo_path <- function(models, data, tree, order = NULL,
                       cor_fun = ape::corPagel, parallel = NULL) {
  cor_fun <- match.fun(cor_fun)
  # Check if all models have the same number of nodes
  var_names <- lapply(models, colnames)
  if (length(models) > 1 &
      (stats::var(lengths(models)) != 0 |
       any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0))) {
    stop('All causal models need to include the same variables.', call. = FALSE)
  }
  if (is.null(names(models))) {
    names(models) <- LETTERS[1:length(models)]
  }
  if ('tbl_df' %in% class(data)) {
    data <- as.data.frame(data)
  }
  if (is.null(order)) {
    order <- find_consensus_order(models)
  }
  formulas <- purrr::map(models, find_formulas, order)
  formulas <- purrr::map(formulas,
                         ~purrr::map(.x, ~{attr(., ".Environment") <- NULL; .}))
  f_list <- unique(unlist(formulas))
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(min(c(parallel::detectCores() - 1,
                                      length(f_list))),
                                parallel)
    parallel::clusterExport(cl, list('gls2', 'data', 'tree', 'cor_fun'),
                            environment())
    on.exit(parallel::stopCluster(cl))
  } else {
    cl <- NULL
  }
  dsep_models <- pbapply::pblapply(f_list, function(x) {
    gls2(x, data = data, tree = tree, cor_fun = cor_fun)
  }, cl = cl)
  dsep_models <- purrr::map(formulas, ~dsep_models[match(.x, f_list)])

  d_sep <- purrr::map2(formulas, dsep_models,
                      ~dplyr::data_frame(d_sep = as.character(.x),
                                         p = purrr::map_dbl(.y, get_p),
                                         corStruct = purrr::map_dbl(.y, ~get_corStruct(.)[[1]]),
                                         model = .y))

  out <- list(d_sep = d_sep, models = models, data = data, tree = tree,
              cor_fun = cor_fun)
  class(out) <- 'phylopath'
  return(out)
}

#' @export
summary.phylopath <- function(object, ...) {
  phylopath <- object
  k <- sapply(phylopath$d_sep, nrow)
  q <- sapply(phylopath$models, function(m) nrow(m) + sum(m))
  C <- sapply(phylopath$d_sep, function(x) C_stat(x$p))
  p <- C_p(C, k)
  IC <- CICc(C, q, nrow(phylopath$data))

  d <- data.frame(model = names(phylopath$models), k = k, q = q, C = C, p = p,
                  CICc = IC, stringsAsFactors = FALSE)
  d <- d[order(d$CICc), ]
  d$delta_CICc <- d$CICc - d$CICc[1]
  d$l <- l(d$delta_CICc)
  d$w <- w(d$l)
  class(d) <- c('phylopath_summary', 'data.frame')
  return(d)
}

#' Extract and estimate the best supported model from a phylogenetic path
#' analysis.
#'
#' @param phylopath An object of class \code{phylopath}.
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   best_model <- best(p)
#'   # Print the best model to see coefficients, se and ci:
#'   best_model
#'   # Plot to show the weighted graph:
#'   plot(best_model)
#'
best <- function(phylopath) {
  b <- summary(phylopath)[1, 'model']
  best_model <- phylopath$models[[b]]
  est_DAG(best_model, phylopath$data, phylopath$cor_fun, phylopath$tree)
}

#' Extract and estimate an arbitrary model from a phylogenetic path analysis.
#'
#' @param phylopath An object of class \code{phylopath}.
#' @param choice A character string of the name of the model to be chosen, or
#'   the index in \code{models}.
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   my_model <- choice(p, "B")
#'   # Print the best model to see coefficients, se and ci:
#'   my_model
#'   # Plot to show the weighted graph:
#'   plot(my_model)
#'
choice <- function(phylopath, choice) {
  est_DAG(phylopath$models[[choice]], phylopath$data, phylopath$cor_fun,
          phylopath$tree)
}

#' Extract and average the best supported models from a phylogenetic path
#' analysis.
#'
#' @param phylopath An object of class \code{phylopath}.
#' @param cut_off The CICc cut-off used to select the best models. Use
#'   \code{Inf} to average over all models. Use the \code{best} function to
#'   only use the top model.
#' @inheritParams average_DAGs
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS),
#'                      C = DAG(LS ~ BM, NL ~ LS + BM, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   summary(p)
#'   # Models A and C have close to equal support, so we may decide to take
#'   # their average.
#'
#'   avg_model <- average(p)
#'   # Print the average model to see coefficients, se and ci:
#'   avg_model
#'
#'   \dontrun{
#'   # Plot to show the weighted graph:
#'   plot(avg_model)
#'   # Note that coefficents that only occur in one of the models become much
#'   # smaller when we use full averaging:
#'   coef_plot(avg_model)
#'   coef_plot(average(p, method = 'full'))
#'   }
#'
average <- function(phylopath, cut_off = 2, method = 'conditional', ...) {
  d <- summary(phylopath)
  b <- d[d$delta_CICc < cut_off, ]
  best_models <- lapply(phylopath$models[b$model], est_DAG, phylopath$data,
                        phylopath$cor_fun, phylopath$tree)
  average <- average_DAGs(best_models, b$w, method, ...)

  class(average$coef) <- c('matrix', 'DAG')
  return(average)
}
