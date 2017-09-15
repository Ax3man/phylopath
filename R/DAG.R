#' Directed acyclic graphs (DAGs)
#'
#' This function is a simple wrapper around the function from the \code{ggm}
#' package with the same name. The only differences are that the \code{order}
#' argument defaults to \code{TRUE} and that it adds a \code{DAG} class for
#' easy plotting. Typically, one would use \code{\link{build_model_set}} to
#' create models for use with the \code{phylopath} package.
#'
#' Supply a formulas for the model as arguments. Formulas should be of the
#' form `parent ~ child` and describe each path in your model. Multiple
#' children of a single parent can be combined into a single formula:
#' `parent ~ child1 + child2`. Finally, an isolate (unconnected variable) can
#' be included as being connected to itself: `isolate ~ isolate`.
#'
#' @param order logical, defaulting to `TRUE` If `TRUE` the nodes of the DAG
#'   are permuted according to the topological order. If `FALSE` the nodes are
#'   in the order they first appear in the model formulae (from left to right).
#'   For use in the `phylopath` package, this should always be kept to `TRUE`,
#'   but the argument is available to avoid potential problems with masking the
#'   function from other packages.
#'
#' @inheritParams ggm::DAG
#' @return An object of classes \code{matrix} and \code{DAG}
#' @export
#'
#' @examples
#'   # Use formula notation to create DAGs:
#'   plot(DAG(A~B, B~C))
#'   # Use + to easily add multiple parents to a node:
#'   plot(DAG(A~B+C))
#'   # Add a node as it's own parent to create an isolate:
#'   plot(DAG(A~B+C, D~D))
DAG <- function(..., order = TRUE) {
  d <- ggm::DAG(..., order = order)
  class(d) <- c(class(d), 'DAG')
  d
}

#' Build a model set.
#'
#' This is a convenience function to quickly and clearly define a set of causal
#' models. Supply a list of formulas for each model, using either `list()`,
#' or `c()`. Formulas should be of the form `parent ~ child` and describe each
#' path in your model. Multiple children of a single parent can be combined
#' into a single formula: `parent ~ child1 + child2`.
#'
#' @section Note This function uses `ggm::DAG`
#'
#' @param ... Named arugments, which each are a lists of formulas defining the
#'   paths of a causal model.
#' @param .common A list of formulas that contain causal paths that are common
#'   to each model.
#'
#' @return A list of models, each of class \code{matrix} and \code{DAG}.
#' @export
#'
#' @examples
#' (m <- build_model_set(
#'   A = c(a~b, b~c),
#'   B = c(b~a, c~b),
#'   .common = c(d~a)))
#' plot_model_set(m)
build_model_set <- function(..., .common) {
  model_list <- list(...)
  # Get all unique variables
  vars <- unique(unlist(lapply(unlist(model_list), all.vars)))
  # And guarantee their inclusion as isolates if necessary
  vars_formulas <- lapply(vars, function(x) stats::as.formula(paste(x, '~', x)))
  .common <- c(.common, vars_formulas, recursive = TRUE)
  # Add isolates and common paths to all models
  model_list <- lapply(model_list, function(x) c(x, .common))
  # Build the models with DAG
  lapply(model_list, function(x) do.call(DAG, x))
}

#' Add standardized path coefficients to a DAG.
#'
#' @param DAG A directed acyclic graph, typically created with \code{DAG}.
#' @inheritParams phylo_path
#'
#' @return An object of class \code{fitted_DAG}.
#'
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   plot(d)
#'   d_fitted <- est_DAG(d, rhino, ape::corBrownian, rhino_tree)
#'   plot(d_fitted)
est_DAG <- function(DAG, data, cor_fun, tree) {
  cor_fun <- match.fun(cor_fun)
  r <- rownames(data)
  data <- dplyr::mutate_if(data, is.numeric, scale)
  rownames(data) <- r
  d <- Map(function(x, y, n) {
    if (all(y == 0)) {
      return(cbind(y, y, y, y))
    }
    f <- stats::formula(paste(x, paste(n[y == 1], collapse = '+'), sep = '~'))
    m <- gls2(f, data = data, cor_fun = cor_fun, tree = tree)
    if (!is.null(m$error)) {
      stop(paste('Fitting the following model:\n   ', Reduce(paste, deparse(f)),
                 '\nproduced this error:\n   ', m$error),
           call. = FALSE)
    }
    m <- m$result
    Coef <- se <- lower <- upper <- y
    Coef[Coef != 0]   <- get_est(m)
    se[se != 0]       <- get_se(m)
    lower[lower != 0] <- tryCatch(get_lower(m),
                                  error = function(e) NA)
    upper[upper != 0] <- tryCatch(get_upper(m),
                                  error = function(e) NA)
    return(cbind(coef = Coef, se = se, lower = lower, upper = upper))
  }, colnames(DAG), as.data.frame(DAG), MoreArgs = list(n = rownames(DAG)))
  if (any(sapply(d, function(x) any(is.na(x))))) {
    warnings("NA's have been generated, most likely some confidence intervals could not be estimated.")
  }
  coefs  <- sapply(d, `[`, 1:nrow(DAG), 1)
  ses    <- sapply(d, `[`, 1:nrow(DAG), 2)
  lowers <- sapply(d, `[`, 1:nrow(DAG), 3)
  uppers <- sapply(d, `[`, 1:nrow(DAG), 4)
  rownames(coefs) <- rownames(ses) <- rownames(lowers) <- rownames(uppers) <-
    rownames(DAG)
  res <- list(coef = coefs, se = ses, lower = lowers, upper = uppers)
  class(res) <- 'fitted_DAG'
  return(res)
}

#' Add standardized path coefficients to a binary DAG.
#'
#' @param DAG A directed acyclic graph, typically created with \code{DAG}.
#' @inheritParams phylo_path_binary
#'
#' @return An object of class \code{binary_fitted_DAG}.
#'
#' @export
est_DAG_binary <- function(DAG, data, tree) {
  d <- Map(function(x, y, n) {
    if (all(y == 0)) {
      return(cbind(y, y, y, y))
    }
    f <- stats::formula(paste(x, paste(n[y == 1], collapse = '+'), sep = '~'))
    m <- purrr::safely(ape::binaryPGLMM)(f, data = data, phy = tree)
    if (!is.null(m$error)) {
      stop(paste('Fitting the following model:\n   ', Reduce(paste, deparse(f)),
                 '\nproduced this error:\n   ', m$error),
           call. = FALSE)
    }
    m <- m$result
    Coef <- se <- y
    Coef[Coef != 0]   <- get_est_binary(m)
    se[se != 0]       <- get_se_binary(m)
    return(cbind(coef = Coef, se = se))
  }, colnames(DAG), as.data.frame(DAG), MoreArgs = list(n = rownames(DAG)))
  if (any(sapply(d, function(x) any(is.na(x))))) {
    warnings("NA's have been generated, most likely some confidence intervals could not be estimated.")
  }
  coefs  <- sapply(d, `[`, 1:nrow(DAG), 1)
  ses    <- sapply(d, `[`, 1:nrow(DAG), 2)
  rownames(coefs) <- rownames(ses) <- rownames(DAG)
  res <- list(coef = coefs, se = ses)
  class(res) <- c('binary_fitted_DAG', 'fitted_DAG')
  return(res)
}

#' Perform model averaging on a list of DAGs.
#'
#' @param fitted_DAGs A list of \code{fitted_DAG} objects containing
#'   coefficients and standard errors, usually obtained by using \code{est_DAG}
#'   on several DAGs.
#' @param weights A vector of associated model weights.
#' @param method Either \code{"full"} or \code{"conditional"}. The methods
#'   differ in how they deal with averaging a path coefficient where the path is
#'   absent in some of the models. The full method sets the coefficient (and the
#'   variance) for the missing paths to zero, meaning paths that are missing in
#'   some models will shrink towards zero. The conditional method only averages
#'   over models where the path appears, making it more sensitive to small
#'   effects. Following von Hardenberg & Gonzalez-Voyer 2013, conditional
#'   averaging is set as the default. Also see \link[MuMIn]{model.avg}.
#' @param ... Additional arguments passed to \link[MuMIn]{par.avg}.
#'
#'   For details on the error calculations, see \link[MuMIn]{par.avg}.
#'
#' @return An object of class \code{fitted_DAG}, including standard errors and
#'   confidence intervals.
#' @export
#'
#' @examples
#'   # Normally, I would advocate the use of the phylo_path and average
#'   # functions, but this code shows how to average any set of models. Note
#'   # that not many checks are implemented, so you may want to be careful and
#'   # make sure the DAGs make sense and contain the same variables!
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   fit_cand <- lapply(candidates, est_DAG, rhino, ape::corPagel, rhino_tree)
#'   ave_cand <- average_DAGs(fit_cand)
#'   coef_plot(ave_cand)
average_DAGs <- function(fitted_DAGs, weights = rep(1, length(coef)),
                         method = 'conditional', ...) {
  if (!(method %in% c('full', 'conditional'))) {
    stop('method has to be either "full" or "conditional".', call. = FALSE)
  }
  ord <- rownames(fitted_DAGs[[1]]$coef)
  fitted_DAGs <- lapply(fitted_DAGs, function(l) {
    lapply(l, function(m) m[ord, ord]) } )

  coef      <- lapply(fitted_DAGs, `[[`, 'coef')
  std_error <- lapply(fitted_DAGs, `[[`, 'se')

  rel_weights <- weights / sum(weights)
  coef      <- simplify2array(coef)
  std_error <- simplify2array(std_error)
  if (method == 'conditional') {
    coef[coef == 0]          <- NA
    std_error[std_error == 0] <- NA
  }
  a_coef <- apply(coef, 1:2, stats::weighted.mean, w = rel_weights,
                  na.rm = TRUE)
  a_coef[is.nan(a_coef)] <- 0

  if (!is.null(std_error)) {
    coef_list      <- purrr::array_branch(coef, 1:2)
    std_error_list <- purrr::array_branch(std_error, 1:2)
    r <- purrr::map2(coef_list, std_error_list,
                     function(.x, .y, ...) MuMIn::par.avg(.x, .y, rel_weights, ...),
                     ...)
    r <- purrr::map(r, function(x) { x[is.nan(x)] <- 0; x } )
    a_std_error <- matrix(purrr::map_dbl(r, "SE"), nrow = nrow(coef))
    lower       <- matrix(purrr::map_dbl(r, "Lower CI"), nrow = nrow(coef))
    upper       <- matrix(purrr::map_dbl(r, "Upper CI"), nrow = nrow(coef))
    dimnames(a_std_error) <- dimnames(lower) <- dimnames(upper) <-
      dimnames(a_coef)
    res <- list(coef = a_coef, se = a_std_error, lower = lower, upper = upper)
  } else {
    res <- list(coef = a_coef)
  }
  class(res) <- 'fitted_DAG'
  return(res)
}