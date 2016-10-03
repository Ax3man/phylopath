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

#' Add standardized path coefficients to a DAG.
#'
#' @param DAG A directed acyclic graph, typically created with \code{DAG}.
#' @inheritParams phylo_path
#'
#' @return An object of class \code{fitted_DAG}.
#'
#' @export
est_DAG <- function(DAG, data, cor_fun, tree) {
  r <- rownames(data)
  data <- dplyr::mutate_if(data, is.numeric, scale)
  rownames(data) <- r
  d <- Map(function(x, y, n) {
    if (all(y == 0)) {
      return(cbind(y, y, y, y))
    }
    f <- stats::formula(paste(x, paste(n[y == 1], collapse = '+'), sep = '~'))
    m <- gls2(f, data = data, cor_fun = cor_fun, tree = tree)
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
  res <- list(coef = coefs, se = ses, lower = lowers, upper = uppers)
  class(res) <- 'fitted_DAG'
  return(res)
}

#' Perform model averaging on a list of DAGs.
#'
#' @param coef A list of estimated DAGs, usually created with \code{est_DAG}.
#' @param std_error Optionally, a list of associated matrices standard errors,
#'   one matrix for each DAG. Can also be obtained from \code{est_DAG}.
#' @param weights A vector of associated model weights.
#' @param method Either \code{"full"} or \code{"conditional"}. The methods
#'   differ in how they deal with averaging a path coefficient where the path is
#'   absent in some of the models. The full method sets the coefficent (and the
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
#'   confidence inverals if \code{std_error} was supplied.
#' @export
average_DAGs <- function(coef, std_error = NULL,
                         weights = rep(1, length(coef)), method = 'conditional',
                         ...) {
  if (!(method %in% c('full', 'conditional'))) {
    stop('method has to be either "full" or "conditional".', call. = FALSE)
  }
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
                     ~ MuMIn::par.avg(.x, .y, rel_weights, ...))
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