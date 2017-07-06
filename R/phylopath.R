#' Compare causal models in a phylogenetic context.
#'
#' @param models A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{DAG}.
#' @param data A \code{data.frame} with data.
#' @param tree A phylogenetic tree of class \code{pylo}.
#' @param cor_fun A function that creates a \code{corStruct} object, typically
#'   one of the \code{cor*} functions from the \code{ape} package, such as
#'   \code{corBrownian}, \code{corPagel} etc.
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
#' @param na.rm Should rows that contain missing values be dropped from the data
#'   as necessary (with a message)?
#'
#' @return A phylopath object, with the following components:
#'  \describe{
#'   \item{d_sep}{for each model a table with separation statements and statistics.}
#'   \item{models}{the DAGs}
#'   \item{data}{the supplied data}
#'   \item{tree}{the supplied tree}
#'   \item{cor_fun}{the employed correlation structure}
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
phylo_path <- function(models, data, tree, cor_fun = ape::corPagel,
                       order = NULL, parallel = NULL, na.rm = TRUE) {
  cor_fun <- match.fun(cor_fun)
  tmp <- check_models_data_tree(models, data, tree, na.rm)
  models <- tmp$models
  data <- tmp$data
  tree <- tmp$tree

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
  dsep_models_runs <- pbapply::pblapply(f_list, function(x) {
    gls2(x, data = data, tree = tree, cor_fun = cor_fun)
  }, cl = cl)
  # Produce appropriate error if needed
  errors <- purrr::map(dsep_models_runs, 'error')
  purrr::map2(errors, f_list,
              ~if(!is.null(.x))
                stop(paste('Fitting the following model:\n   ',
                           Reduce(paste, deparse(f_list[[1]])),
                           '\nproduced this error:\n   ', .x),
                     call. = FALSE))
  # Otherwise collect models and move on.
  dsep_models <- purrr::map(dsep_models_runs, 'result')
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

#' Compare binary causal models in a phylogenetic context.
#'
#' @param models A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{DAG}.
#' @param data A \code{data.frame} with binary data. Variables can be a numeric vector (with all
#'   values either 0 or 1), a factor (with two levels) or character vector (with two unique values).
#' @param tree A phylogenetic tree of class \code{pylo}.
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
#' @param na.rm Should rows that contain missing values be dropped from the data
#'   as necessary (with a message)?
#'
#' @details Estimation of relationships between binary variables is built upon
#'   \code{ape::binaryPGLMM} and stands or falls with the accuracy of that method. Model fitting
#'   using this method is particularly slow and the \code{parallel} argument can therefore be
#'   especially useful.
#'
#' @return A phylopath_binary object, with the following components:
#'  \describe{
#'   \item{d_sep}{for each model a table with separation statements and statistics.}
#'   \item{models}{the DAGs}
#'   \item{data}{the supplied data}
#'   \item{tree}{the supplied tree}
#'   }
#' @export
#'
phylo_path_binary <- function(models, data, tree, order = NULL, parallel = NULL, na.rm = TRUE) {
  tmp <- check_models_data_tree(models, data, tree, na.rm)
  models <- tmp$models
  data <- tmp$data
  tree <- tmp$tree

  # Check if all variables are binary
  data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)
  data[sapply(data, is.factor)] <- lapply(data[sapply(data, is.factor)],
                                          function(x) as.numeric(x) - 1)
  if (!all(as.matrix(data) %in% c(0, 1))) {
    stop("All variables need to binary. Supply them either as vectors of 0 and 1, two level factors,
         or character vectors with two unique values.", call. = FALSE)
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
    parallel::clusterExport(cl, list('data', 'tree'), environment())
    on.exit(parallel::stopCluster(cl))
  } else {
    cl <- NULL
  }
  dsep_models_runs <- pbapply::pblapply(f_list, function(x) {
    purrr::safely(ape::binaryPGLMM)(x, data = data, phy = tree)
  }, cl = cl)
  # Produce appropriate error if needed
  errors <- purrr::map(dsep_models_runs, 'error')
  purrr::map2(errors, f_list,
              ~if(!is.null(.x))
                stop(paste('Fitting the following model:\n   ',
                           Reduce(paste, deparse(f_list[[1]])),
                           '\nproduced this error:\n   ', .x),
                     call. = FALSE))
  # Otherwise collect models and move on.
  dsep_models <- purrr::map(dsep_models_runs, 'result')
  dsep_models <- purrr::map(formulas, ~dsep_models[match(.x, f_list)])

  d_sep <- purrr::map2(formulas, dsep_models,
                       ~dplyr::data_frame(d_sep = as.character(.x),
                                          p = purrr::map_dbl(.y, get_p_binary),
                                          s2 = purrr::map_dbl(.y, 's2'),
                                          model = .y))

  out <- list(d_sep = d_sep, models = models, data = data, tree = tree)
  class(out) <- c('phylopath_binary', 'phylopath')
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
  UseMethod('best')
}

#' @export
best.phylopath <- function(phylopath) {
  b <- summary(phylopath)[1, 'model']
  best_model <- phylopath$models[[b]]
  est_DAG(best_model, phylopath$data, phylopath$cor_fun, phylopath$tree)
}

#' @export
best.phylopath_binary <- function(phylopath) {
  b <- summary(phylopath)[1, 'model']
  best_model <- phylopath$models[[b]]
  est_DAG_binary(best_model, phylopath$data, phylopath$tree)
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
  UseMethod('choice')
}

#' @export
choice.phylopath <- function(phylopath, choice) {
  est_DAG(phylopath$models[[choice]], phylopath$data, phylopath$cor_fun,
          phylopath$tree)
}

#' @export
choice.phylopath_binary <- function(phylopath, choice) {
  est_DAG_binary(phylopath$models[[choice]], phylopath$data, phylopath$tree)
}


#' Extract and average the best supported models from a phylogenetic path
#' analysis.
#'
#' @param phylopath An object of class \code{phylopath} or
#'   \code{phylopath_binary}.
#' @param cut_off The CICc cut-off used to select the best models. Use
#'   \code{Inf} to average over all models. Use the \code{best} function to
#'   only use the top model, or \code{choice} to select any single model.
#' @inheritParams average_DAGs
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(
#'     A = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL),
#'     B = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#'   )
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   summary(p)
#'
#'   # Models A and B have similar support, so we may decide to take
#'   # their average.
#'
#'   avg_model <- average(p)
#'   # Print the average model to see coefficients, se and ci:
#'   avg_model
#'
#'   \dontrun{
#'   # Plot to show the weighted graph:
#'   plot(avg_model)
#'
#'   # One can see that an averaged model is not necessarily a DAG itself.
#'   # This model actually has a path in two directions.
#'
#'   # Note that coefficients that only occur in one of the models become much
#'   # smaller when we use full averaging:
#'
#'   coef_plot(avg_model)
#'   coef_plot(average(p, method = 'full'))
#'   }
#'
average <- function(phylopath, cut_off = 2, method = 'conditional', ...) {
  UseMethod('average')
}

#' @export
average.phylopath <- function(phylopath, cut_off = 2, method = 'conditional', ...) {
  d <- summary(phylopath)
  b <- d[d$delta_CICc < cut_off, ]
  best_models <- lapply(phylopath$models[b$model], est_DAG, phylopath$data,
                        phylopath$cor_fun, phylopath$tree)
  average <- average_DAGs(best_models, b$w, method, ...)
  class(average$coef) <- c('matrix', 'DAG')
  return(average)
}

#' @export
average.phylopath_binary <- function(phylopath, cut_off = 2, method = 'conditional', ...) {
  d <- summary(phylopath)
  b <- d[d$delta_CICc < cut_off, ]
  best_models <- lapply(phylopath$models[b$model], est_DAG_binary, phylopath$data, phylopath$tree)
  average <- average_DAGs(best_models, b$w, method, ...)
  class(average$coef) <- c('matrix', 'DAG')
  return(average)
}
