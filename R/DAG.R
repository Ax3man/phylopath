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

#' Estimate path coefficients for a DAG.
#'
#' This function will estimate the path coefficients for a DAG, and return them
#' with uncertainty.
#'
#' @details
#' Continuous variables are standardized before fitting, that is, centered and
#' divided by their standard deviation, while binary variables are left as they
#' are. In the two usual cases this means the coefficients can be compared with
#' each other directly:
#'
#' * When every variable is continuous, each coefficient is the change in the
#'   outcome, in standard deviations, for an increase of one standard deviation in
#'   the predictor. These are standardized regression coefficients, as reported by
#'   von Hardenberg & Gonzalez-Voyer.
#' * When every variable is binary, each coefficient is the log odds ratio for the
#'   outcome between the two levels of the predictor. A level to level contrast is
#'   a natural unit that every binary variable has, so these are comparable too.
#'
#' When a model contains binary *and* continuous variables, its coefficients can no
#' longer be compared with one another, for two separate reasons. Paths into a
#' binary variable are log odds ratios while paths into a continuous variable are
#' standardized regression coefficients, which are different units entirely. And a
#' binary predictor contributes a contrast between its two levels, which for a
#' variable where a proportion `p` of species has one of the levels amounts to
#' `1 / sqrt(p * (1 - p))` standard deviations. That is always more than two, so
#' the coefficients of binary predictors are inflated relative to those of
#' continuous predictors in the same model.
#'
#' @param DAG A directed acyclic graph, typically created with \code{DAG}.
#' @param boot The number of bootstrap replicates used to estimate confidence intervals.
#' @inheritParams phylo_path
#'
#' @return An object of class \code{fitted_DAG}, which is a list with the
#'   following components:
#'   \describe{
#'    \item{coef}{a matrix of path coefficients, with predictors in the rows and
#'      responses in the columns. Absent paths are 0.}
#'    \item{se}{a matrix of standard errors of the coefficients.}
#'    \item{lower, upper}{matrices with the bounds of the bootstrapped
#'      confidence intervals, only present if `boot` was larger than 0.}
#'    \item{binary}{a named logical vector marking which variables were modelled
#'      as binary.}
#'   }
#'
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   plot(d)
#'   d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
#'   plot(d_fitted)
est_DAG <- function(DAG, data, tree, model = 'lambda', method = 'logistic_MPLE',
                    boot = 0, ...) {
  stopifnot(inherits(DAG, 'DAG'))
  dots <- check_dots(list(...))
  # check and convert binary variables
  vars <- rownames(DAG)
  data[vars] <- as_binary_factors(data[vars])
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
      rlang::abort(c(
        'A phylogenetic regression could not be fitted.',
        x = paste('Model:', Reduce(paste, deparse(f))),
        x = paste('Error:', m$error)
      ))
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
  # Record which variables were modelled as binary, since that determines what
  # the coefficients of paths into them mean.
  res$binary <- vapply(data[rownames(DAG)], is.factor, logical(1))
  class(res) <- 'fitted_DAG'
  warn_if_mixed(res)
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
#' @seealso [est_DAG()] for what the coefficients mean.
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
  if (length(list(...)) != 0) {
    rlang::abort(
      c('The `...` argument of `average_DAGs()` is deprecated.',
        x = paste0('Given: ', paste(names(list(...)), collapse = ', '), '.'))
    )
  }
  avg_method <- match.arg(avg_method, choices = c("full", "conditional"))
  ord <- rownames(fitted_DAGs[[1]]$coef)
  binary <- combine_binary(fitted_DAGs, ord)

  fitted_DAGs <- lapply(fitted_DAGs, function(l) {
    mats <- intersect(c('coef', 'se', 'lower', 'upper'), names(l))
    l[mats] <- lapply(l[mats], function(m) m[ord, ord])
    l
  } )

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
  res$binary <- binary
  class(res) <- 'fitted_DAG'
  warn_if_mixed(res)
  return(res)
}

#' Extract the paths of a fitted causal model.
#'
#' Returns the estimated paths of a `fitted_DAG` as a long `data.frame`, one row
#' per path, which is a convenient starting point for custom tables and figures.
#' Paths that are absent from the causal model are not included.
#'
#' @param x An object of class `fitted_DAG`.
#' @param row.names,optional,... Ignored, present for consistency with the
#'   generic.
#' @param order_by Either `"default"`, to order the paths by their position in
#'   the causal model, `"causal"`, to follow the paths downstream, or
#'   `"strength"`, to order them by the size of their coefficients. Ordering by
#'   strength is not available for models that contain both binary and continuous
#'   variables, since their coefficients are not comparable, see [est_DAG()].
#'   Ordering by cause is not available for cyclical models, which [average()] can
#'   produce when the direction of a path is not resolved.
#'
#' @return A `data.frame` with columns `from`, `to`, `coef` and `se`, as well as
#'   `lower` and `upper` if the model was fitted with `boot` larger than 0.
#' @export
#'
#' @examples
#'   d <- DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS)
#'   d_fitted <- est_DAG(d, rhino, rhino_tree, 'lambda')
#'   as.data.frame(d_fitted)
#'   as.data.frame(d_fitted, order_by = 'strength')
as.data.frame.fitted_DAG <- function(x, row.names = NULL, optional = FALSE, ...,
                                     order_by = 'default') {
  order_by <- match.arg(order_by, c('default', 'causal', 'strength'), FALSE)
  if (order_by == 'strength' && isTRUE(is_mixed(x))) {
    rlang::abort(
      c('`order_by = "strength"` orders the paths by the size of their coefficients, which is not meaningful for a causal model that contains both binary and continuous variables.',
        x = 'Those coefficients are not on comparable scales, see `?est_DAG`.',
        i = 'Use `order_by = "causal"` or `"default"` instead.')
    )
  }
  coef <- x$coef

  ind <- which(coef != 0, arr.ind = TRUE)
  out <- data.frame(
    from = rownames(coef)[ind[, 'row']],
    to   = colnames(coef)[ind[, 'col']],
    coef = coef[ind],
    se   = x$se[ind]
  )
  if (!is.null(x$lower)) {
    out$lower <- x$lower[ind]
    out$upper <- x$upper[ind]
  }
  vars <- switch(
    order_by,
    default  = rownames(coef),
    causal   = causal_order(coef),
    strength = NULL
  )
  out <- if (is.null(vars)) out[order(out$coef), ]
         else out[order(match(out$from, vars), match(out$to, vars)), ]
  rownames(out) <- NULL
  out
}
