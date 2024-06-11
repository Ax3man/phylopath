#' Directed acyclic graphs (DAGs)
#'
#' This function is a simple wrapper around the function from the `ggm`
#' package with the same name. The only differences are that the `order`
#' argument defaults to `TRUE` and that it adds a `DAG` class for
#' easy plotting. Typically, one would use [define_model_set()] to
#' create models for use with the `phylopath` package.
#'
#' Supply a formulas for the model as arguments. Formulas should be of the
#' form `child ~ parent`` and describe each path in your model. Multiple
#' children of a single parent can be combined into a single formula:
#' `child ~ parent1 + parent2`. Finally, an isolate (unconnected variable) can
#' be included as being connected to itself: `isolate ~ isolate`.
#'
#' @param order logical, defaulting to `TRUE`. If `TRUE` the nodes of the DAG
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

#' Define a model set.
#'
#' This is a convenience function to quickly and clearly define a set of causal
#' models. Supply a list of formulas for each model, using either `c()`. Formulas
#' should be of the form `child ~ parent` and describe each path in your model.
#' Multiple children of a single parent can be combined into a single formula:
#' `child ~ parent1 + parent2`.
#'
#' This function uses [ggm::DAG()].
#'
#' @param ... Named arguments, which each are a lists of formulas defining the
#'   paths of a causal model.
#' @param .common A list of formulas that contain causal paths that are common
#'   to each model.
#'
#' @return A list of models, each of class `matrix` and `DAG`.
#' @export
#'
#' @examples
#' (m <- define_model_set(
#'   A = c(a~b, b~c),
#'   B = c(b~a, c~b),
#'   .common = c(d~a)))
#' plot_model_set(m)
define_model_set <- function(..., .common = NULL) {
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
#' @param boot The number of bootstrap replicates used to estimate confidence intervals.
#' @inheritParams phylo_path
#'
#' @return An object of class \code{fitted_DAG}.
#'
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   plot(d)
#'   d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
#'   plot(d_fitted)
est_DAG <- function(DAG, data, tree, model, method, boot = 0, ...) {
  stopifnot(inherits(DAG, 'DAG'))
  dots <- list(...)
  # scale the continous variables
  r <- rownames(data)
  data[sapply(data, is.numeric)] <- lapply(data[sapply(data, is.numeric)], scale)
  rownames(data) <- r
  d <- Map(function(x, y, n) {
    if (all(y == 0)) {
      return(cbind(y, y, y, y))
    }
    f <- stats::formula(paste(x, paste(n[y == 1], collapse = '+'), sep = '~'))
    m <- phylo_g_lm(f, data, tree, model, method, boot, dots)
    if (!is.null(m$error)) {
      stop(paste('Fitting the following model:\n   ', Reduce(paste, deparse(f)),
                 '\nproduced this error:\n   ', m$error),
           call. = FALSE)
    }
    m <- m$result
    Coef <- se <- lower <- upper <- y
    Coef[Coef != 0]   <- get_est(m)
    se[se != 0]       <- get_se(m)
    lower[lower != 0] <- get_lower(m)
    upper[upper != 0] <- get_upper(m)
    return(cbind(coef = Coef, se = se, lower = lower, upper = upper))
  }, colnames(DAG), as.data.frame(DAG), MoreArgs = list(n = rownames(DAG)))
  coefs  <- sapply(d, `[`, 1:nrow(DAG), 1)
  ses    <- sapply(d, `[`, 1:nrow(DAG), 2)
  lowers <- sapply(d, `[`, 1:nrow(DAG), 3)
  uppers <- sapply(d, `[`, 1:nrow(DAG), 4)
  rownames(coefs) <- rownames(ses) <- rownames(lowers) <- rownames(uppers) <-
    rownames(DAG)
  if (boot > 0) {
    res <- list(coef = coefs, se = ses, lower = lowers, upper = uppers)
  } else {
    res <- list(coef = coefs, se = ses)
  }
  class(res) <- 'fitted_DAG'
  return(res)
}

#' Perform model averaging on a list of DAGs.
#'
#' @param fitted_DAGs A list of `fitted_DAG` objects containing
#'   coefficients and standard errors, usually obtained by using [est_DAG()]
#'   on several DAGs.
#' @param weights A vector of associated model weights.
#' @param avg_method Either `"full"` or `"conditional"`. The methods
#'   differ in how they deal with averaging a path coefficient where the path is
#'   absent in some of the models. The full method sets the coefficient (and the
#'   variance) for the missing paths to zero, meaning paths that are missing in
#'   some models will shrink towards zero. The conditional method only averages
#'   over models where the path appears, making it more sensitive to small
#'   effects. Following von Hardenberg & Gonzalez-Voyer 2013, conditional
#'   averaging is set as the default.
#' @param ... Use of the ellipses is deprecated.
#'
#'
#' @return An object of class `fitted_DAG`, including standard errors and
#'   confidence intervals.
#' @export
#'
#' @examples
#'   # Normally, I would advocate the use of the phylo_path and average
#'   # functions, but this code shows how to average any set of models. Note
#'   # that not many checks are implemented, so you may want to be careful and
#'   # make sure the DAGs make sense and contain the same variables!
#'   candidates <- define_model_set(
#'     A = NL ~ BM,
#'     B = NL ~ LS,
#'     .common = c(LS ~ BM, DD ~ NL)
#'   )
#'   fit_cand <- lapply(candidates, est_DAG, rhino, rhino_tree,
#'                      model = 'lambda', method = 'logistic_MPLE')
#'   ave_cand <- average_DAGs(fit_cand)
#'   coef_plot(ave_cand)
average_DAGs <- function(fitted_DAGs, weights = rep(1, length(coef)),
                         avg_method = 'conditional', ...) {
  if (length(list(...)) != 0) stop('Use of ... is deprecated.')
  avg_method <- match.arg(avg_method, choices = c("full", "conditional"))
  ord <- rownames(fitted_DAGs[[1]]$coef)
  fitted_DAGs <- lapply(fitted_DAGs, function(l) {
    lapply(l, function(m) m[ord, ord]) } )

  coef      <- lapply(fitted_DAGs, `[[`, 'coef')
  std_error <- lapply(fitted_DAGs, `[[`, 'se')

  rel_weights <- weights / sum(weights)
  coef      <- simplify2array(coef)
  std_error <- simplify2array(std_error)
  if (avg_method == 'conditional') {
    coef[coef == 0]           <- NA
    std_error[std_error == 0] <- NA
  }
  a_coef <- apply(coef, 1:2, stats::weighted.mean, w = rel_weights, na.rm = TRUE)
  a_coef[is.nan(a_coef)] <- 0

  if (!is.null(std_error)) {
    coef_list      <- purrr::array_branch(coef, 1:2)
    std_error_list <- purrr::array_branch(std_error, 1:2)
    r <- purrr::map2(
      coef_list,
      std_error_list,
      function(.x, .y, ...) par_avg(.x, .y, rel_weights, ...),
      ...
    )
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